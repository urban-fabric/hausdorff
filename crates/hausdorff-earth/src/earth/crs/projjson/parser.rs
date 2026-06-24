/// A minimal parsed JSON value — enough to read PROJJSON, no external crate.
/// `Null`/`Bool` complete the JSON model (PROJJSON has e.g. `deprecated` booleans)
/// even though this reader never inspects their contents.
pub(super) enum Json {
    Null,
    #[allow(dead_code)]
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    /// Parse a JSON document, or `None` if malformed.
    pub(super) fn parse(s: &str) -> Option<Json> {
        let mut p = JsonParser {
            b: s.as_bytes(),
            i: 0,
        };
        let v = p.value()?;
        p.ws();
        // Trailing garbage is tolerated (PROJJSON documents are single values).
        Some(v)
    }

    /// The value for `key` in an object.
    pub(super) fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(m) => m.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub(super) fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    pub(super) fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(n) => Some(*n),
            _ => None,
        }
    }

    pub(super) fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(a) => Some(a),
            _ => None,
        }
    }
}

/// A byte-cursor recursive-descent JSON parser.
struct JsonParser<'a> {
    b: &'a [u8],
    i: usize,
}

impl JsonParser<'_> {
    fn ws(&mut self) {
        while matches!(self.b.get(self.i), Some(b' ' | b'\n' | b'\t' | b'\r')) {
            self.i += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.i).copied()
    }

    fn value(&mut self) -> Option<Json> {
        self.ws();
        match self.peek()? {
            b'{' => self.object(),
            b'[' => self.array(),
            b'"' => Some(Json::Str(self.string()?)),
            b't' => self.lit("true", Json::Bool(true)),
            b'f' => self.lit("false", Json::Bool(false)),
            b'n' => self.lit("null", Json::Null),
            _ => self.number(),
        }
    }

    fn lit(&mut self, word: &str, val: Json) -> Option<Json> {
        if self.b[self.i..].starts_with(word.as_bytes()) {
            self.i += word.len();
            Some(val)
        } else {
            None
        }
    }

    fn object(&mut self) -> Option<Json> {
        self.i += 1; // '{'
        let mut m = Vec::new();
        self.ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Some(Json::Obj(m));
        }
        loop {
            self.ws();
            let key = self.string()?;
            self.ws();
            if self.peek()? != b':' {
                return None;
            }
            self.i += 1;
            m.push((key, self.value()?));
            self.ws();
            match self.peek()? {
                b',' => self.i += 1,
                b'}' => {
                    self.i += 1;
                    return Some(Json::Obj(m));
                }
                _ => return None,
            }
        }
    }

    fn array(&mut self) -> Option<Json> {
        self.i += 1; // '['
        let mut a = Vec::new();
        self.ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Some(Json::Arr(a));
        }
        loop {
            a.push(self.value()?);
            self.ws();
            match self.peek()? {
                b',' => self.i += 1,
                b']' => {
                    self.i += 1;
                    return Some(Json::Arr(a));
                }
                _ => return None,
            }
        }
    }

    fn string(&mut self) -> Option<String> {
        if self.peek()? != b'"' {
            return None;
        }
        self.i += 1;
        let mut s = String::new();
        loop {
            let c = self.peek()?;
            self.i += 1;
            match c {
                b'"' => return Some(s),
                b'\\' => {
                    let e = self.peek()?;
                    self.i += 1;
                    match e {
                        b'"' => s.push('"'),
                        b'\\' => s.push('\\'),
                        b'/' => s.push('/'),
                        b'n' => s.push('\n'),
                        b't' => s.push('\t'),
                        b'r' => s.push('\r'),
                        b'b' => s.push('\u{8}'),
                        b'f' => s.push('\u{c}'),
                        b'u' => {
                            let hex = std::str::from_utf8(self.b.get(self.i..self.i + 4)?).ok()?;
                            let cp = u32::from_str_radix(hex, 16).ok()?;
                            self.i += 4;
                            s.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                        }
                        _ => return None,
                    }
                }
                // ASCII, or a raw UTF-8 byte we never inspect (PROJJSON is ASCII-escaped).
                _ => s.push(c as char),
            }
        }
    }

    fn number(&mut self) -> Option<Json> {
        let start = self.i;
        while matches!(self.peek(), Some(c)
            if c.is_ascii_digit() || matches!(c, b'-' | b'+' | b'.' | b'e' | b'E'))
        {
            self.i += 1;
        }
        std::str::from_utf8(&self.b[start..self.i])
            .ok()?
            .parse()
            .ok()
            .map(Json::Num)
    }
}
