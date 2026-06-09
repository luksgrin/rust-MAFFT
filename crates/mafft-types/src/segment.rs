/// An aligned segment region.
///
/// Replaces the C `Segment` struct which uses raw pointer-based pairing.
/// Here, the paired segment is referenced by an optional index into the
/// segment collection.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentSegment {
    pub start: i32,
    pub end: i32,
    pub center: i32,
    pub score: f64,
    pub skip_forward: i32,
    pub skip_backward: i32,
    /// Index of the paired segment, if any.
    pub pair: Option<usize>,
    pub number: i32,
}

impl Default for AlignmentSegment {
    fn default() -> Self {
        Self {
            start: 0,
            end: 0,
            center: 0,
            score: 0.0,
            skip_forward: 0,
            skip_backward: 0,
            pair: None,
            number: 0,
        }
    }
}

/// A pair of segments from two groups being aligned.
///
/// Replaces the C `Segments` struct.
#[derive(Debug, Clone)]
pub struct SegmentPair {
    pub group1: AlignmentSegment,
    pub group2: AlignmentSegment,
    pub number1: i32,
    pub number2: i32,
}
