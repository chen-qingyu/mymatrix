use mymatrix::{Matrix, MatrixError, Vector};
use pyinrs::Fraction;
use rstest::{fixture, rstest};

struct Fixture {
    mat_0x0: Matrix,
    mat_1x1: Matrix,
    mat_3x3: Matrix,
}

#[fixture]
fn setup() -> Fixture {
    Fixture {
        mat_0x0: Matrix::new(),
        mat_1x1: Matrix::from([[2]]),
        mat_3x3: Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
    }
}

#[rstest]
fn basics(setup: Fixture) {
    assert_eq!(setup.mat_0x0.row_size(), 0);
    assert_eq!(setup.mat_0x0.col_size(), 0);
    assert!(setup.mat_0x0.is_empty());

    assert_eq!(setup.mat_1x1.row_size(), 1);
    assert_eq!(setup.mat_1x1.col_size(), 1);
    assert!(!setup.mat_1x1.is_empty());

    assert_eq!(setup.mat_3x3.row_size(), 3);
    assert_eq!(setup.mat_3x3.col_size(), 3);
    assert!(!setup.mat_3x3.is_empty());

    assert_eq!(Matrix::from([[1, 2, 3]]).row_size(), 1);
    assert_eq!(Matrix::from([[1, 2, 3]]).col_size(), 3);
    assert!(!Matrix::from([[1, 2, 3]]).is_empty());
}

#[rstest]
fn is_square(setup: Fixture) {
    assert!(setup.mat_0x0.is_square());
    assert!(setup.mat_1x1.is_square());
    assert!(setup.mat_3x3.is_square());
    assert!(!Matrix::from([[1, 2, 3]]).is_square());
    assert!(!Matrix::from([[1], [2]]).is_square());
}

#[rstest]
fn from_vec_vector() {
    let vecs = vec![Vector::from([1, 2, 3]), Vector::from([4, 5, 6])];
    let m = Matrix::from(vecs);
    assert_eq!(m, Matrix::from([[1, 2, 3], [4, 5, 6]]));
}

#[rstest]
fn compare(setup: Fixture) {
    assert!(setup.mat_3x3 == Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]));
    assert!(setup.mat_3x3 != Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]));
}

#[rstest]
fn access(mut setup: Fixture) {
    for r in 0..setup.mat_3x3.row_size() {
        for c in 0..setup.mat_3x3.col_size() {
            assert_eq!(setup.mat_3x3[r][c], ((r * 3 + c + 1) as i32).into());
        }
    }

    setup.mat_3x3[0][0] = 0.into();
    assert_eq!(setup.mat_3x3[0][0], 0.into());
}

#[rstest]
fn is_symmetric(setup: Fixture) {
    assert!(setup.mat_0x0.is_symmetric());
    assert!(setup.mat_1x1.is_symmetric());
    assert!(!setup.mat_3x3.is_symmetric());

    assert!(Matrix::identity(3).is_symmetric());
}

#[rstest]
fn is_upper(setup: Fixture) {
    assert!(setup.mat_0x0.is_upper());
    assert!(setup.mat_1x1.is_upper());
    assert!(!setup.mat_3x3.is_upper());

    assert!(Matrix::identity(3).is_upper());
}

#[rstest]
fn is_lower(setup: Fixture) {
    assert!(setup.mat_0x0.is_lower());
    assert!(setup.mat_1x1.is_lower());
    assert!(!setup.mat_3x3.is_lower());

    assert!(Matrix::identity(3).is_lower());
}

#[rstest]
fn is_diagonal(setup: Fixture) {
    assert!(setup.mat_0x0.is_diagonal());
    assert!(setup.mat_1x1.is_diagonal());
    assert!(!setup.mat_3x3.is_diagonal());

    assert!(Matrix::identity(3).is_diagonal());
}

#[rstest]
fn trace(setup: Fixture) {
    assert_eq!(setup.mat_0x0.trace(), 0.into());
    assert_eq!(setup.mat_1x1.trace(), 2.into());
    assert_eq!(setup.mat_3x3.trace(), 15.into());

    assert_eq!(Matrix::identity(3).trace(), 3.into());
}

#[rstest]
fn diag() {
    // non-square: min(row, col) elements
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).diag(), Vector::from([1, 5]));
    assert_eq!(Matrix::identity(3).diag(), Vector::ones(3));
    assert_eq!(Matrix::from([[7]]).diag(), Vector::from([7]));
    assert_eq!(Matrix::new().diag(), Vector::new());

    // round trip
    let v = Vector::from([2, 3, 4]);
    let d = Matrix::from_diagonal(&v);
    assert_eq!(d, Matrix::from([[2, 0, 0], [0, 3, 0], [0, 0, 4]]));
    assert_eq!(d.diag(), v);
}

#[rstest]
fn transpose(setup: Fixture) {
    assert_eq!(setup.mat_0x0.transpose(), setup.mat_0x0);
    assert_eq!(setup.mat_1x1.transpose(), setup.mat_1x1);
    assert_eq!(setup.mat_3x3.transpose().transpose(), setup.mat_3x3);

    assert_eq!(Matrix::zeros(2, 3).transpose(), Matrix::zeros(3, 2));
}

#[rstest]
fn col() {
    let m = Matrix::from([[1, 2, 3], [4, 5, 6]]);
    assert_eq!(m.col(0), Vector::from([1, 4]));
    assert_eq!(m.col(1), Vector::from([2, 5]));
    assert_eq!(m.col(2), Vector::from([3, 6]));

    // consistency with row access via transpose
    assert_eq!(m.col(1), m.transpose()[1]);
}

#[rstest]
fn rotate_left(setup: Fixture) {
    assert_eq!(setup.mat_0x0.rotate_left(), Matrix::new());
    assert_eq!(setup.mat_1x1.rotate_left(), setup.mat_1x1);
    assert_eq!(setup.mat_3x3.rotate_left(), Matrix::from([[3, 6, 9], [2, 5, 8], [1, 4, 7]]));

    // 4 rotations = identity
    assert_eq!(setup.mat_3x3.rotate_left().rotate_left().rotate_left().rotate_left(), setup.mat_3x3);

    // non-square
    assert_eq!(
        Matrix::from([[1, 2, 3, 4], [5, 6, 7, 8]]).rotate_left(),
        Matrix::from([[4, 8], [3, 7], [2, 6], [1, 5]])
    );
}

#[rstest]
fn rotate_right(setup: Fixture) {
    assert_eq!(setup.mat_0x0.rotate_right(), Matrix::new());
    assert_eq!(setup.mat_1x1.rotate_right(), setup.mat_1x1);
    assert_eq!(setup.mat_3x3.rotate_right(), Matrix::from([[7, 4, 1], [8, 5, 2], [9, 6, 3]]));

    // 4 rotations = identity
    assert_eq!(setup.mat_3x3.rotate_right().rotate_right().rotate_right().rotate_right(), setup.mat_3x3);

    // rotate_left = rotate_right * 3
    assert_eq!(setup.mat_3x3.rotate_right().rotate_right().rotate_right(), setup.mat_3x3.rotate_left());

    // non-square
    assert_eq!(
        Matrix::from([[1, 2, 3, 4], [5, 6, 7, 8]]).rotate_right(),
        Matrix::from([[5, 1], [6, 2], [7, 3], [8, 4]])
    );
}

#[rstest]
fn kron() {
    let a = Matrix::from([[1, 2], [3, 4]]);
    let b = Matrix::from([[0, 5], [6, 7]]);
    assert_eq!(a.kron(&b), Matrix::from([[0, 5, 0, 10], [6, 7, 12, 14], [0, 15, 0, 20], [18, 21, 24, 28],]));
    // A ⊗ I: block diagonal
    assert_eq!(a.kron(&Matrix::identity(2)), Matrix::from([[1, 0, 2, 0], [0, 1, 0, 2], [3, 0, 4, 0], [0, 3, 0, 4]]));
    // mixed sizes: (2x1) ⊗ (1x2)
    assert_eq!(Matrix::from([[1], [2]]).kron(&Matrix::from([[3, 4]])), Matrix::from([[3, 4], [6, 8]]));
}

#[rstest]
fn hadamard() {
    let a = Matrix::from([[1, 2], [3, 4]]);
    let b = Matrix::from([[5, 6], [7, 8]]);
    assert_eq!(a.hadamard(&b), Matrix::from([[5, 12], [21, 32]]));
    assert_eq!(a.hadamard(&Matrix::ones(2, 2)), a);
}

#[rstest]
fn row_echelon_form(setup: Fixture) {
    assert_eq!(setup.mat_0x0.row_echelon_form(), Matrix::new());
    assert_eq!(setup.mat_1x1.row_echelon_form(), Matrix::from([[2]]));
    assert_eq!(setup.mat_3x3.row_echelon_form(), Matrix::from([[1, 2, 3], [0, -3, -6], [0, 0, 0]]));

    assert_eq!(Matrix::ones(2, 2).row_echelon_form(), Matrix::from([[1, 1], [0, 0]]));
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).row_echelon_form(), Matrix::from([[1, 2, 3], [0, -3, -6]]));
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).row_echelon_form(), Matrix::from([[1, 2], [0, -2], [0, 0]]));
}

#[rstest]
fn row_canonical_form(setup: Fixture) {
    assert_eq!(setup.mat_0x0.row_canonical_form(), Matrix::new());
    assert_eq!(setup.mat_1x1.row_canonical_form(), Matrix::from([[1]]));
    assert_eq!(setup.mat_3x3.row_canonical_form(), Matrix::from([[1, 0, -1], [0, 1, 2], [0, 0, 0]]));

    assert_eq!(Matrix::ones(2, 2).row_canonical_form(), Matrix::from([[1, 1], [0, 0]]));
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).row_canonical_form(), Matrix::from([[1, 0, -1], [0, 1, 2]]));
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).row_canonical_form(), Matrix::from([[1, 0], [0, 1], [0, 0]]));
}

#[rstest]
fn det(setup: Fixture) {
    assert_eq!(setup.mat_0x0.det(), 1.into());
    assert_eq!(setup.mat_1x1.det(), 2.into());
    assert_eq!(setup.mat_3x3.det(), 0.into());

    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).det(), 27.into());
}

#[rstest]
fn submatrix(setup: Fixture) {
    assert_eq!(setup.mat_3x3.submatrix(0, 0), Matrix::from([[5, 6], [8, 9]]));
    assert_eq!(setup.mat_3x3.submatrix(1, 1), Matrix::from([[1, 3], [7, 9]]));
    assert_eq!(setup.mat_3x3.submatrix(0, 2), Matrix::from([[4, 5], [7, 8]]));
}

#[rstest]
#[should_panic(expected = "Error: Index out of range.")]
fn bad_submatrix_empty() {
    Matrix::new().submatrix(0, 0);
}

#[rstest]
fn minor(setup: Fixture) {
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).minor(), Matrix::from([[4, 3], [2, 1]]));
    assert_eq!(setup.mat_3x3.minor(), Matrix::from([[-3, -6, -3], [-6, -12, -6], [-3, -6, -3]]));
}

#[rstest]
fn cofactor(setup: Fixture) {
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).cofactor(), Matrix::from([[4, -3], [-2, 1]]));
    assert_eq!(setup.mat_3x3.cofactor(), Matrix::from([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]]));
}

#[rstest]
fn adj(setup: Fixture) {
    assert_eq!(setup.mat_0x0.adj(), Matrix::new());
    assert_eq!(setup.mat_1x1.adj(), Matrix::from([[1]]));
    assert_eq!(setup.mat_3x3.adj(), Matrix::from([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]]));
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).adj(), Matrix::from([[4, -2], [-3, 1]]));

    // properties
    let m = &setup.mat_3x3;
    assert_eq!(Matrix::identity(3).adj(), Matrix::identity(3)); // I.adj = I
    assert_eq!(Matrix::zeros(3, 3).adj(), Matrix::zeros(3, 3)); // O.adj = O (n>1)
    assert_eq!(Matrix::zeros(1, 1).adj(), Matrix::identity(1)); // O.adj = I (n=1)
    assert_eq!((2 * m.clone()).adj(), 4 * m.adj()); // (cA).adj = c^(n-1)*A.adj
    assert_eq!(m.transpose().adj(), m.adj().transpose()); // (A^T).adj = (A.adj)^T
    assert_eq!(m.adj().det(), m.det() * m.det()); // |A.adj| = |A|^(n-1)
    assert_eq!(m * m.adj(), m.det() * Matrix::identity(3)); // A*A.adj = |A|*I
    assert_eq!((m * m).adj(), m.adj() * m.adj()); // (A^k).adj = (A.adj)^k
    assert_eq!(m.adj().adj(), m.det() * m.clone()); // A.adj.adj = |A|^(n-2)*A
}

#[rstest]
fn inv(setup: Fixture) {
    assert_eq!(setup.mat_0x0.inv(), Ok(Matrix::new()));
    assert_eq!(setup.mat_1x1.inv(), Ok(Matrix::from([[Fraction::from((1, 2))]])));
    assert_eq!(setup.mat_3x3.inv(), Err(MatrixError::Singular));

    assert_eq!(
        Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).inv(),
        Ok(Matrix::from([
            [Fraction::from((-16, 9)), Fraction::from((8, 9)), Fraction::from((-1, 9))],
            [Fraction::from((14, 9)), Fraction::from((-7, 9)), Fraction::from((2, 9))],
            [Fraction::from((-1, 9)), Fraction::from((2, 9)), Fraction::from((-1, 9))],
        ]))
    );
}

#[rstest]
fn rank(setup: Fixture) {
    assert_eq!(setup.mat_0x0.rank(), 0);
    assert_eq!(setup.mat_1x1.rank(), 1);
    assert_eq!(setup.mat_3x3.rank(), 2);

    assert_eq!(Matrix::ones(2, 2).rank(), 1);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).rank(), 3);
}

#[rstest]
fn subspaces(setup: Fixture) {
    // rank 2, nullity 1
    let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    assert_eq!(m.row_space(), Matrix::from([[1, 2, 3], [0, -3, -6]]));
    assert_eq!(m.col_space(), Matrix::from([[1, 4, 7], [2, 5, 8]]));
    assert_eq!(m.null_space(), Matrix::from([[1, -2, 1]]));

    // full-rank square: trivial null space, bases are identity
    assert_eq!(Matrix::identity(3).row_space(), Matrix::identity(3));
    assert_eq!(Matrix::identity(3).col_space(), Matrix::identity(3));
    assert_eq!(Matrix::identity(3).null_space(), Matrix::new());

    // rank 1: ones(2,2)
    let o = Matrix::ones(2, 2);
    assert_eq!(o.row_space(), Matrix::from([[1, 1]]));
    assert_eq!(o.col_space(), Matrix::from([[1, 1]]));
    assert_eq!(o.null_space(), Matrix::from([[-1, 1]]));

    // rectangular: 3x2 full column rank
    let a = Matrix::from([[1, 2], [3, 4], [5, 6]]);
    assert_eq!(a.row_space(), Matrix::from([[1, 2], [0, -2]]));
    assert_eq!(a.col_space(), Matrix::from([[1, 3, 5], [2, 4, 6]]));
    assert_eq!(a.null_space(), Matrix::new());

    // empty
    assert_eq!(setup.mat_0x0.null_space(), Matrix::new());
    assert_eq!(setup.mat_0x0.row_space(), Matrix::new());
    assert_eq!(setup.mat_0x0.col_space(), Matrix::new());
}

#[rstest]
fn characteristic_polynomial() {
    // 1x1: λ - a
    assert_eq!(Matrix::from([[3]]).characteristic_polynomial(), vec![1.into(), (-3).into()]);
    // 2x2: λ² - tr·λ + det
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).characteristic_polynomial(), vec![1.into(), (-5).into(), (-2).into()]);
    // diagonal 3x3: (λ-2)(λ-3)(λ-4)
    assert_eq!(
        Matrix::from([[2, 0, 0], [0, 3, 0], [0, 0, 4]]).characteristic_polynomial(),
        vec![1.into(), (-9).into(), 26.into(), (-24).into()]
    );

    // verify p(λ) == det(λI - A) at λ = 2
    let a = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]);
    let p = a.characteristic_polynomial();
    let lam = 2;
    let det = Matrix::from([[lam - 1, -2, -3], [-4, lam - 5, -6], [-7, -8, lam]]).det();
    let mut val = Fraction::new();
    for &c in &p {
        val = val * Fraction::from(lam) + c;
    }
    assert_eq!(val, det);
}

#[rstest]
fn pseudo_inverse() {
    // invertible -> equals the inverse
    let a = Matrix::from([[1, 2], [3, 4]]);
    assert_eq!(a.pseudo_inverse(), a.inv().unwrap());

    // singular rank-1
    let s = Matrix::from([[1, 2], [2, 4]]);
    assert_eq!(
        s.pseudo_inverse(),
        Matrix::from([[Fraction::from((1, 25)), Fraction::from((2, 25))], [Fraction::from((2, 25)), Fraction::from((4, 25))],])
    );
    // Moore-Penrose condition: A A⁺ A = A
    assert_eq!(&s * &s.pseudo_inverse() * &s, s);

    // non-square full column rank: A⁺ A = I
    let r = Matrix::from([[1, 2], [3, 4], [5, 6]]);
    let rp = r.pseudo_inverse();
    assert_eq!(&r * &rp * &r, r);
    assert_eq!(&rp * &r, Matrix::identity(2));

    // zero matrix -> zero matrix
    assert_eq!(Matrix::zeros(2, 3).pseudo_inverse(), Matrix::zeros(3, 2));
}

#[rstest]
fn lu_decomposition(setup: Fixture) {
    assert_eq!(
        Matrix::from([[2, 3, 1], [4, 7, 1], [6, 7, 3]]).lu_decomposition(),
        Ok((Matrix::from([[1, 0, 0], [2, 1, 0], [3, -2, 1]]), Matrix::from([[2, 3, 1], [0, 1, -1], [0, 0, -2]])))
    );

    // singular matrix -> Err(Singular)
    assert_eq!(setup.mat_3x3.lu_decomposition(), Err(MatrixError::Singular));
}

#[rstest]
fn cholesky() {
    // A = [[4, 2], [2, 3]] is SPD
    // Expected: L = [[1, 0], [1/2, 1]], D = [4, 2]
    let a = Matrix::from([[4, 2], [2, 3]]);
    let (l, d) = a.cholesky().unwrap();
    assert_eq!(l, Matrix::from([[Fraction::from(1), Fraction::from(0)], [Fraction::from((1, 2)), Fraction::from(1)]]));
    assert_eq!(d, Vector::from([4, 2]));

    // verify A = L * D * L^T
    let n = a.row_size();
    let mut d_mat = Matrix::zeros(n, n);
    for i in 0..n {
        d_mat[i][i] = d[i];
    }
    assert_eq!(&l * &d_mat * &l.transpose(), a);

    // 3x3 example: A = [[5, 4, 2], [4, 5, 2], [2, 2, 3]]
    let a3 = Matrix::from([[5, 4, 2], [4, 5, 2], [2, 2, 3]]);
    let (l3, d3) = a3.cholesky().unwrap();
    let mut d3_mat = Matrix::zeros(3, 3);
    for i in 0..3 {
        d3_mat[i][i] = d3[i];
    }
    assert_eq!(&l3 * &d3_mat * &l3.transpose(), a3);

    // non-symmetric -> NotPositiveDefinite
    let non_sym = Matrix::from([[1, 2], [3, 4]]);
    assert_eq!(non_sym.cholesky(), Err(MatrixError::NotPositiveDefinite));

    // non-positive-definite -> NotPositiveDefinite
    let non_pd = Matrix::from([[-1, 0], [0, -1]]);
    assert_eq!(non_pd.cholesky(), Err(MatrixError::NotPositiveDefinite));

    // identity
    let (l_id, d_id) = Matrix::identity(3).cholesky().unwrap();
    assert_eq!(l_id, Matrix::identity(3));
    assert_eq!(d_id, Vector::ones(3));
}

#[rstest]
fn pow() {
    let a = Matrix::from([[1, 2], [3, 4]]);
    assert_eq!(a.pow(0), Matrix::identity(2));
    assert_eq!(a.pow(1), a);
    // A^2 = [[7, 10], [15, 22]]
    assert_eq!(a.pow(2), Matrix::from([[7, 10], [15, 22]]));
    // A^3 = A^2 * A = [[37, 54], [81, 118]]
    assert_eq!(a.pow(3), Matrix::from([[37, 54], [81, 118]]));
}

#[rstest]
fn solve() {
    // 2x + 3y = 7, 4x + 5y = 13 => x=2, y=1
    let a = Matrix::from([[2, 3], [4, 5]]);
    let b = Vector::from([7, 13]);
    assert_eq!(a.solve(&b), Ok(Vector::from([2, 1])));

    // consistent but singular (infinitely many solutions)
    let singular = Matrix::from([[1, 2], [2, 4]]);
    let b2 = Vector::from([3, 6]);
    assert_eq!(singular.solve(&b2), Err(MatrixError::Singular));

    // inconsistent system (no solution)
    let m = Matrix::from([[1, 1], [1, 1]]);
    assert_eq!(m.solve(&Vector::from([2, 3])), Err(MatrixError::Inconsistent));
}

#[rstest]
fn general_solution() {
    // unique solution: null space is empty
    let a = Matrix::from([[2, 3], [4, 5]]);
    let (xp, null) = a.general_solution(&Vector::from([7, 13])).unwrap();
    assert_eq!(xp, Vector::from([2, 1]));
    assert_eq!(null, Matrix::new());

    // underdetermined but consistent: x + 2y = 3 (rank 1)
    let s = Matrix::from([[1, 2], [2, 4]]);
    let (xp, null) = s.general_solution(&Vector::from([3, 6])).unwrap();
    assert_eq!(xp, Vector::from([3, 0]));
    assert_eq!(null, Matrix::from([[-2, 1]]));
    assert_eq!(&s * &xp, Vector::from([3, 6]));
    assert_eq!(&s * &null[0], Vector::zeros(2));

    // x + y = 2
    let m = Matrix::from([[1, 1], [1, 1]]);
    let (xp, null) = m.general_solution(&Vector::from([2, 2])).unwrap();
    assert_eq!(xp, Vector::from([2, 0]));
    assert_eq!(null, Matrix::from([[-1, 1]]));
    assert_eq!(&m * &xp, Vector::from([2, 2]));
    assert_eq!(&m * &null[0], Vector::zeros(2));

    // inconsistent -> Err(Inconsistent)
    assert_eq!(m.general_solution(&Vector::from([2, 3])), Err(MatrixError::Inconsistent));
}

#[rstest]
fn split() {
    let matrix = Matrix::from([[1, 2], [3, 4], [5, 6]]);

    assert_eq!(matrix.split_row(1).0, Matrix::from([[1, 2]]));
    assert_eq!(matrix.split_row(1).1, Matrix::from([[3, 4], [5, 6]]));

    assert_eq!(matrix.split_col(1).0, Matrix::from([[1], [3], [5]]));
    assert_eq!(matrix.split_col(1).1, Matrix::from([[2], [4], [6]]));
}

#[rstest]
fn expand() {
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]).expand_row(Matrix::zeros(2, 2)),
        &Matrix::from([[1, 2], [3, 4], [0, 0], [0, 0]])
    );
    assert_eq!(
        Matrix::from([[1, 2, 3, 4, 5]]).expand_row(Matrix::from([[6, 7, 8, 9, 10]])),
        &Matrix::from([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
    );

    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]).expand_col(Matrix::zeros(2, 2)),
        &Matrix::from([[1, 2, 0, 0], [3, 4, 0, 0]])
    );
    assert_eq!(
        Matrix::from([[1, 2, 3, 4, 5]]).expand_col(Matrix::from([[6, 7, 8, 9, 10]])),
        &Matrix::from([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]])
    );
}

#[rstest]
fn elementary_row_operations(mut setup: Fixture) {
    assert_eq!(setup.mat_3x3.e_row_swap(0, 1), &Matrix::from([[4, 5, 6], [1, 2, 3], [7, 8, 9]]));
    assert_eq!(setup.mat_3x3.e_scalar_multiplication(1, 2.into()), &Matrix::from([[4, 5, 6], [2, 4, 6], [7, 8, 9]]));
    assert_eq!(setup.mat_3x3.e_row_sum(0, 1, (-1).into()), &Matrix::from([[2, 1, 0], [2, 4, 6], [7, 8, 9]]));
}

#[rstest]
fn add() {
    assert_eq!(Matrix::create(2, 3, 1.into()) + Matrix::create(2, 3, 2.into()), Matrix::create(2, 3, 3.into()));
    assert_eq!(Matrix::from([[1, 2], [3, 4]]) + Matrix::from([[-1, -1], [-1, -1]]), Matrix::from([[0, 1], [2, 3]]));
}

#[rstest]
fn sub() {
    assert_eq!(Matrix::create(2, 3, 1.into()) - Matrix::create(2, 3, 2.into()), Matrix::create(2, 3, (-1).into()));
    assert_eq!(Matrix::from([[1, 2], [3, 4]]) - Matrix::from([[-1, -1], [-1, -1]]), Matrix::from([[2, 3], [4, 5]]));
}

#[rstest]
fn mul() {
    assert_eq!(Matrix::create(2, 2, 1.into()) * Matrix::create(2, 2, 2.into()), Matrix::create(2, 2, 4.into()));
    assert_eq!(Matrix::create(1, 3, 1.into()) * Matrix::create(3, 1, 1.into()), Matrix::create(1, 1, 3.into()));

    assert_eq!(Matrix::create(2, 3, 1.into()) * Fraction::from(2), Matrix::create(2, 3, 2.into()));
    assert_eq!(Matrix::from([[1, 2], [3, 4]]) * Fraction::from(3), Matrix::from([[3, 6], [9, 12]]));
}

#[rstest]
fn scalar_div() {
    assert_eq!(Matrix::create(2, 2, 2.into()) / Fraction::from(2), Matrix::create(2, 2, 1.into()));
    assert_eq!(Matrix::from([[3, 6], [9, 12]]) / Fraction::from(3), Matrix::from([[1, 2], [3, 4]]));
}

#[rstest]
fn mul_vector() {
    let m = Matrix::from([[1, 2], [3, 4], [5, 6]]);
    let v = Vector::from([1, 2]);
    assert_eq!(&m * &v, Vector::from([5, 11, 17]));

    // identity
    assert_eq!(&Matrix::identity(3) * &Vector::from([7, 8, 9]), Vector::from([7, 8, 9]));

    // zero matrix
    assert_eq!(&Matrix::zeros(2, 3) * &Vector::from([1, 2, 3]), Vector::zeros(2));
}

#[rstest]
fn format(setup: Fixture) {
    assert_eq!(format!("{}", setup.mat_0x0), "[\n]");
    assert_eq!(format!("{}", setup.mat_1x1), "[\n2\n]");
    assert_eq!(format!("{}", setup.mat_3x3), "[\n1 2 3\n4 5 6\n7 8 9\n]");

    assert_eq!(
        format!(
            "{}",
            Matrix::from([[Fraction::from((-11, 6)), Fraction::from((5, 6))], [Fraction::from((5, 3)), Fraction::from((-2, 3))]])
        ),
        "[
-11/6   5/6
  5/3  -2/3
]"
    );
}

#[rstest]
fn from_iterator() {
    let m: Matrix = [vec![1, 2], vec![3, 4]].into_iter().collect();
    assert_eq!(m, Matrix::from([[1, 2], [3, 4]]));

    let rows: Matrix = [Vector::from([1, 2]), Vector::from([3, 4])].into_iter().collect();
    assert_eq!(rows, Matrix::from([[1, 2], [3, 4]]));

    let fracs: Matrix = [vec![Fraction::from(1), Fraction::from(2)]].into_iter().collect();
    assert_eq!(fracs, Matrix::from([[1, 2]]));
}

#[rstest]
fn iter_methods() {
    let m = Matrix::from([[1, 2], [3, 4]]);
    assert_eq!(m.iter().cloned().collect::<Matrix>(), m);

    // IntoIterator for &Matrix
    let mut nrows = 0;
    for row in &m {
        assert_eq!(row.size(), 2);
        nrows += 1;
    }
    assert_eq!(nrows, 2);

    // IntoIterator for &mut Matrix
    let mut w = Matrix::from([[1, 2], [3, 4]]);
    for row in &mut w {
        row[0] = 0.into();
    }
    assert_eq!(w, Matrix::from([[0, 2], [0, 4]]));
}
