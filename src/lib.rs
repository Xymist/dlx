use std::collections::HashMap;

type Index = usize;

#[derive(Clone, Debug, Default, PartialEq)]
// Represents a 1 in the sparse matrix, or a header entry. Contains links to nearest
// neighbours plus some header data if appropriate.
struct MtxElem {
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

/// An assembled matrix, for solving a specific problem
struct Solver {
    pub subset_rows: HashMap<usize, usize>,
    pub elements: Vec<MtxElem>,
    pub bottoms: Vec<Index>,

    // The number of columns still participating in the algorithm. If this
    // falls to zero, a solution has been found.
    pub active_columns: usize,
}

impl Solver {
    // For reasons of efficiency, we choose the column containing
    // the smallest number of ones (which in this sparse matrix is
    // the smallest column) at each iteration of the algorithm.
    // This minimises the number of levels in each branch of the tree.
    fn smallest_column(&self) -> Index {
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

    pub fn solve(&mut self) {}

    fn subsolve(&mut self) {}

    fn cover(&mut self, col_idx: Index) {
        // Find the horizontal neighbours of the chosen column
        let r = self.elements[col_idx].right;
        let l = self.elements[col_idx].left;

        // Remove this header by attaching its right to its left
        // and it's left to its right.
        self.elements[r].left = self.elements[col_idx].left;
        self.elements[l].right = self.elements[col_idx].right;

        // Take the next element of this column, 'downwards'
        let mut nxt = self.elements[col_idx].down;

        // Do the same thing by joining the top and bottom for every other one
        // of every row which has a 1 in that column.
        while nxt != col_idx {
            let mut rgt = self.elements[nxt].right;
            while rgt != nxt {
                // Find the header for this column and reduce its count
                // of ones by one
                let header = self.elements[rgt].col;
                self.elements[header].size -= 1;

                let col_up = self.elements[rgt].up;
                self.elements[col_up].down = self.elements[rgt].down;
                let col_down = self.elements[rgt].down;
                self.elements[col_down].up = self.elements[rgt].up;

                rgt = self.elements[rgt].right;
            }
            nxt = self.elements[nxt].down;
        }
    }

    fn uncover(&mut self, col_idx: Index) {
        // Iterate backwards, restoring the rows previously removed.
        // Starts with the 'top' of the header and works upwards,
        // as the cover operation starts with the header and works down.
        let mut nxt = self.elements[col_idx].up;

        while nxt != col_idx {
            // As with the vertical traversal, this is inverted; it
            // travels 'left' around the loop rather than 'right'
            let mut lft = self.elements[nxt].left;
            while lft != nxt {
                // Increase the ones-count of the header by one
                let col_hdr = self.elements[lft].col;
                self.elements[col_hdr].size += 1;

                let col_up = self.elements[lft].up;
                self.elements[col_up].down = lft;

                let col_down = self.elements[lft].down;
                self.elements[col_down].up = lft;

                lft = self.elements[lft].left;
            }

            nxt = self.elements[nxt].up;
        }

        // Restore the column header
        let r = self.elements[col_idx].right;
        self.elements[r].left = col_idx;

        let l = self.elements[col_idx].left;
        self.elements[l].right = col_idx;
    }
}

impl From<Mtx> for Solver {
    fn from(m: Mtx) -> Self {
        Solver {
            subset_rows: m.subset_rows,
            elements: m.elements,
            active_columns: m.bottoms.len(),
            bottoms: m.bottoms,
        }
    }
}

#[derive(Clone, Debug, Default)]
struct Mtx {
    // Matches each row in the matrix with the ID of its subset,
    // for later analysis
    pub subset_rows: HashMap<usize, usize>,
    // The matrix itself
    pub elements: Vec<MtxElem>,
    // The 'bottom' elements of each column; i.e. the populated element
    // in that column with the highest row index
    pub bottoms: Vec<Index>,
}

impl Mtx {
    pub fn size(&self) -> usize {
        self.elements.len()
    }

    pub fn last_row(&self) -> Index {
        match self.elements.get(self.elements.len() - 1) {
            Some(i) => i.row,
            None => 0,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        let matrix = Mtx {
            subset_rows: HashMap::new(),
            elements: Vec::with_capacity(capacity),
            bottoms: Vec::with_capacity(capacity),
        };

        matrix.push_headers(capacity)
    }

    fn push_headers(mut self, column_count: usize) -> Self {
        // Push the first header element
        self.elements.push(MtxElem {
            col: 0,
            right: 1,
            left: column_count - 1,
            ..Default::default()
        });

        // Push the middle N headers, linked left and right
        for idx in 1..(column_count - 1) {
            self.elements.push(MtxElem {
                col: idx,
                right: idx + 1,
                left: idx - 1,
                ..Default::default()
            })
        }

        // Push the final header element and wrap
        self.elements.push(MtxElem {
            col: column_count - 1,
            right: 0,
            left: column_count - 2,
            ..Default::default()
        });

        self.bottoms.extend(0..column_count);

        self
    }

    pub fn push_subsets<T: PartialEq>(
        mut self,
        collection: &[Subset<T>],
        universe: &Universe<T>,
    ) -> Self {
        for subset in collection.iter() {
            self.push_subset_as_row(subset, universe)
        }
        self
    }

    fn push_subset_as_row<T: PartialEq>(&mut self, s: &Subset<T>, u: &Universe<T>) {
        let row = self.last_row() + 1;
        let row0 = self.size();

        // Empty subsets can be ignored, since they don't contribute
        // to the union.
        if s.elements.len() == 0 {
            return;
        };

        for (col, c) in u.iter().enumerate() {
            if s.elements.contains(&c) {
                let idx = self.size();
                let elem = MtxElem {
                    col,
                    row,
                    left: idx - 1,
                    right: idx + 1,
                    up: self.bottoms[col],
                    ..Default::default()
                };

                self.elements.push(elem);

                self.elements[self.bottoms[col]].down = idx;

                self.bottoms[col] = idx;

                self.elements[col].size += 1;
            }
        }
        let rown = self.elements.len() - 1;

        self.elements[row0].left = rown;
        self.elements[rown].right = row0;

        // This allows us to determine from the solution the exact subsets
        // which comprise it, even when they contain the same elements.
        self.subset_rows.entry(row).or_insert(s.id);
    }

    // Connects the 'bottom' of each column to the 'top',
    // turning the matrix from a ring into a torus. This
    // precludes adding more rows
    pub fn finalise(mut self) -> Solver {
        for (i, &b) in self.bottoms.iter().enumerate() {
            self.elements[b].down = i;
            self.elements[i].up = b;
        }
        self.into()
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
}

pub fn weighted_exact_cover<T: PartialEq>(
    u: Universe<T>,
    c: Collection<T>,
) -> Option<Vec<Subset<T>>> {
    let _solver = Mtx::with_capacity(u.len()).push_subsets(&c, &u).finalise();

    return None;
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
