use std::fmt::{self, Write as _};

use super::S2CellId;

impl S2CellId {
    /// A hex token of the id with leading zeros kept and trailing zeros
    /// stripped (e.g. the level-0 face cells are `"1"`, `"3"`, …, `"b"`).
    #[must_use]
    pub fn to_token(self) -> String {
        let s = format!("{:016x}", self.0);
        let trimmed = s.trim_end_matches('0');
        if trimmed.is_empty() {
            "X".to_string()
        } else {
            trimmed.to_string()
        }
    }

    /// Parses a hex token produced by [`S2CellId::to_token`]. Returns `None` if
    /// the token is malformed.
    #[must_use]
    pub fn from_token(s: &str) -> Option<Self> {
        if s.is_empty() || s.len() > 16 {
            return None;
        }
        let mut n = u64::from_str_radix(s, 16).ok()?;
        if s.len() < 16 {
            n <<= 4 * (16 - s.len() as u32);
        }
        Some(Self(n))
    }
}

impl fmt::Display for S2CellId {
    /// Renders as `"<face>/<child positions>"`, e.g. `"1/3210"`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if !self.is_valid() {
            return write!(f, "Invalid:{:x}", self.0);
        }
        write!(f, "{}", self.face())?;
        f.write_char('/')?;
        for level in 1..=self.level() {
            write!(f, "{}", self.child_position(level))?;
        }
        Ok(())
    }
}
