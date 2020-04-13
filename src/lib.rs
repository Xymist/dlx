use std::collections::HashMap;

type Index = usize;

#[derive(Clone, Debug, Default, PartialEq)]
// Represents a 1 in the sparse matrix, or a header entry. Contains links to nearest
// neighbours plus some header data if appropriate.
pub struct MtxElem {
    // The index of the column header, or itself if this is a header
    col: Index,
    // Row index (0 for headers)
    row: Index,
    // The number of elements (1s) in this column (0 for non-headers or empty columns)
    size: Index,
    // The 1 element immediately prior in the row. Circular/wraps.
    left: Index,
    // The 1 element immediately after in the row. Circular/wraps.
    right: Index,
    // The 1 element immediately above in the column. Circular/wraps.
    up: Index,
    // The 1 element immediately below in the column. Circular/wraps.
    down: Index,
}

#[derive(Clone, Debug, Default)]
pub struct Mtx {
    pub subset_rows: HashMap<usize, usize>,
    pub elements: Vec<MtxElem>,
    pub bottoms: Vec<Index>,
    pub active_columns: usize,
}

impl Mtx {
    pub fn last_row(&self) -> Index {
        match self.elements.get(self.elements.len() - 1) {
            Some(i) => i.row,
            None => 0,
        }
    }

    pub fn size(&self) -> usize {
        self.elements.len()
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Mtx {
            subset_rows: HashMap::new(),
            elements: Vec::with_capacity(capacity),
            bottoms: Vec::with_capacity(capacity),
            active_columns: 0,
        }
    }

    pub fn connect_vertical(&mut self) {
        for (i, &b) in self.bottoms.iter().enumerate() {
            self.elements[b].down = i;
            self.elements[i].up = b;
        }
    }

    pub fn smallest_column(&self) -> Index {
        let mut smol_count = Index::max_value();
        let mut smol_idx: Index = 0;
        let mut col_idx: Index = 0;

        while col_idx > 0 {
            let col_size = self.elements[col_idx].size;

            if col_size < smol_count {
                smol_idx = col_idx;
                smol_count = col_size;
            }

            col_idx = self.elements[col_idx].right
        }

        smol_idx
    }
}

type Universe<T> = Vec<T>;
type Collection<T> = Vec<Subset<T>>;

// A strict subset of the universe to be covered, plus some data.
#[derive(Clone, Debug, PartialEq)]
pub struct Subset<T: PartialEq> {
    // Uniquely identifies the Subset
    id: usize,
    // This variation supports weighted subsets,
    // as a decider between multiple exact covers if present
    weight: f64,
    elements: Vec<T>,
}

impl<T: PartialEq> Subset<T> {
    pub fn new(id: usize, weight: f64, elements: Vec<T>) -> Self {
        Subset {
            id,
            weight,
            elements,
        }
    }

    pub fn push_as_row(&self, matrix: &mut Mtx, u: &Universe<T>) {
        let row = matrix.last_row() + 1;
        let row0 = matrix.size();

        // Empty subsets can be ignored, since they don't contribute
        // to the union.
        if self.elements.len() == 0 {
            return;
        };

        for (col, c) in u.iter().enumerate() {
            if self.elements.contains(&c) {
                let idx = matrix.size();
                let elem = MtxElem {
                    col,
                    row,
                    left: idx - 1,
                    right: idx + 1,
                    up: matrix.bottoms[col],
                    ..Default::default()
                };

                matrix.elements.push(elem);

                matrix.elements[matrix.bottoms[col]].down = idx;

                matrix.bottoms[col] = idx;

                matrix.elements[col].size += 1;
            }
        }
        let rown = matrix.elements.len() - 1;

        matrix.elements[row0].left = rown;
        matrix.elements[rown].right = row0;

        // This allows us to determine from the solution the exact subsets
        // which comprise it, even when they contain the same elements.
        matrix.subset_rows.entry(row).or_insert(self.id);
    }
}

pub fn weighted_exact_cover<T: PartialEq>(
    u: Universe<T>,
    c: Collection<T>,
) -> Option<Vec<Subset<T>>> {
    let mut solution = None;

    let mut matrix = Mtx::with_capacity(u.len());
    let column_count = u.len();
    push_headers(&mut matrix, column_count);
    for subset in c.iter() {
        subset.push_as_row(&mut matrix, &u)
    }
    matrix.connect_vertical();

    return solution;
}

fn push_headers(matrix: &mut Mtx, column_count: usize) {
    // Push the first header element
    matrix.elements.push(MtxElem {
        col: 0,
        right: 1,
        left: column_count - 1,
        ..Default::default()
    });

    // Push the middle N headers, linked left and right
    for idx in 1..(column_count - 1) {
        matrix.elements.push(MtxElem {
            col: idx,
            right: idx + 1,
            left: idx - 1,
            ..Default::default()
        })
    }

    // Push the final header element and wrap
    matrix.elements.push(MtxElem {
        col: column_count - 1,
        right: 0,
        left: column_count - 2,
        ..Default::default()
    });

    matrix.bottoms.extend(0..column_count);
    matrix.active_columns = column_count;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_cover() {
        let u: Vec<usize> = vec![1, 2, 3, 4, 5];
        let s: Vec<Subset<usize>> = vec![
            Subset::new(2, 1.0, vec![1, 3]),
            Subset::new(3, 1.0, vec![2, 3]),
            Subset::new(4, 1.0, vec![2, 4]),
        ];

        // No subset here includes '5', so there is no cover, exact or otherwise
        assert_eq!(None, weighted_exact_cover(u, s));
    }

    #[test]
    fn exact_cover() {
        let u: Vec<usize> = vec![1, 2, 3, 4];
        let s: Vec<Subset<usize>> = vec![
            Subset::new(2, 1.0, vec![1, 3]),
            Subset::new(3, 1.0, vec![2, 3]),
            Subset::new(4, 1.0, vec![2, 4]),
        ];

        assert_eq!(
            Some(vec![s[0].clone(), s[2].clone()]),
            weighted_exact_cover(u, s)
        );
    }
}
