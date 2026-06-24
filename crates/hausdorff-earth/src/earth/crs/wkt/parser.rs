/// A parsed WKT node `KEYWORD[args]`. An arg is a quoted string, a number, a
/// bareword (kept as a [`WktArg::Str`]), or a nested node.
pub(super) struct WktNode {
    pub(super) keyword: String,
    args: Vec<WktArg>,
}

enum WktArg {
    Str(String),
    Num(f64),
    Node(WktNode),
}

impl WktNode {
    /// Parse a full WKT string into its root node.
    pub(super) fn parse(s: &str) -> Option<WktNode> {
        let mut p = WktParser {
            b: s.as_bytes(),
            i: 0,
        };
        p.skip_ws();
        let node = p.node()?;
        Some(node)
    }

    /// The first direct child node with this keyword (case-insensitive).
    pub(super) fn child(&self, keyword: &str) -> Option<&WktNode> {
        self.args.iter().find_map(|a| match a {
            WktArg::Node(n) if n.keyword.eq_ignore_ascii_case(keyword) => Some(n),
            _ => None,
        })
    }

    /// All direct child nodes with this keyword (case-insensitive).
    pub(super) fn children<'a>(&'a self, keyword: &'a str) -> impl Iterator<Item = &'a WktNode> {
        self.args.iter().filter_map(move |a| match a {
            WktArg::Node(n) if n.keyword.eq_ignore_ascii_case(keyword) => Some(n),
            _ => None,
        })
    }

    /// The first node with this keyword anywhere in the subtree (this node first).
    pub(super) fn find_deep(&self, keyword: &str) -> Option<&WktNode> {
        if self.keyword.eq_ignore_ascii_case(keyword) {
            return Some(self);
        }
        self.args.iter().find_map(|a| match a {
            WktArg::Node(n) => n.find_deep(keyword),
            _ => None,
        })
    }

    /// The `i`-th positional arg as a string (quoted string or bareword).
    pub(super) fn str_arg(&self, i: usize) -> Option<&str> {
        match self.args.get(i)? {
            WktArg::Str(s) => Some(s),
            _ => None,
        }
    }

    /// The `i`-th positional arg as a number.
    pub(super) fn num_arg(&self, i: usize) -> Option<f64> {
        match self.args.get(i)? {
            WktArg::Num(n) => Some(*n),
            _ => None,
        }
    }

    /// The EPSG code from this node's direct `ID["EPSG",code]` (or `AUTHORITY`).
    pub(super) fn epsg(&self) -> Option<u32> {
        self.children("ID")
            .chain(self.children("AUTHORITY"))
            .find(|id| {
                id.str_arg(0)
                    .is_some_and(|a| a.eq_ignore_ascii_case("EPSG"))
            })
            .and_then(|id| id.num_arg(1))
            .map(|n| n as u32)
    }
}

/// A byte-cursor recursive-descent parser over a WKT string.
struct WktParser<'a> {
    b: &'a [u8],
    i: usize,
}

impl WktParser<'_> {
    fn skip_ws(&mut self) {
        while matches!(self.b.get(self.i), Some(b' ' | b'\n' | b'\t' | b'\r')) {
            self.i += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.b.get(self.i).copied()
    }

    /// A keyword / bareword token: `[A-Za-z0-9_]+`.
    fn word(&mut self) -> Option<String> {
        let start = self.i;
        while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == b'_') {
            self.i += 1;
        }
        if self.i == start {
            return None;
        }
        Some(String::from_utf8_lossy(&self.b[start..self.i]).into_owned())
    }

    /// A node `KEYWORD[args]` — assumes the keyword starts at the cursor.
    fn node(&mut self) -> Option<WktNode> {
        let keyword = self.word()?;
        self.skip_ws();
        if self.peek() != Some(b'[') {
            return None;
        }
        let args = self.args()?;
        Some(WktNode { keyword, args })
    }

    /// The bracketed, comma-separated argument list — assumes the cursor is at `[`.
    fn args(&mut self) -> Option<Vec<WktArg>> {
        self.i += 1; // consume '['
        let mut args = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Some(args);
        }
        loop {
            self.skip_ws();
            args.push(self.arg()?);
            self.skip_ws();
            match self.peek()? {
                b',' => self.i += 1,
                b']' => {
                    self.i += 1;
                    return Some(args);
                }
                _ => return None,
            }
        }
    }

    /// One argument: quoted string, number, nested node, or bareword.
    fn arg(&mut self) -> Option<WktArg> {
        match self.peek()? {
            b'"' => Some(WktArg::Str(self.quoted())),
            c if c == b'-' || c == b'+' || c == b'.' || c.is_ascii_digit() => {
                let start = self.i;
                while matches!(self.peek(), Some(c)
                    if c.is_ascii_digit() || matches!(c, b'.' | b'-' | b'+' | b'e' | b'E'))
                {
                    self.i += 1;
                }
                std::str::from_utf8(&self.b[start..self.i])
                    .ok()?
                    .parse()
                    .ok()
                    .map(WktArg::Num)
            }
            c if c.is_ascii_alphabetic() => {
                let word = self.word()?;
                self.skip_ws();
                if self.peek() == Some(b'[') {
                    Some(WktArg::Node(WktNode {
                        keyword: word,
                        args: self.args()?,
                    }))
                } else {
                    Some(WktArg::Str(word)) // bareword (e.g. `ellipsoidal`, `north`)
                }
            }
            _ => None,
        }
    }

    /// A double-quoted string, with WKT's `""` → `"` escape. Assumes the cursor
    /// is at the opening quote.
    fn quoted(&mut self) -> String {
        self.i += 1; // opening quote
        let mut s = String::new();
        while let Some(c) = self.peek() {
            self.i += 1;
            if c == b'"' {
                if self.peek() == Some(b'"') {
                    self.i += 1;
                    s.push('"');
                } else {
                    break;
                }
            } else {
                s.push(c as char);
            }
        }
        s
    }
}
