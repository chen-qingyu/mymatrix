use mymatrix::Matrix;
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
fn transpose(setup: Fixture) {
    assert_eq!(setup.mat_0x0.transpose(), setup.mat_0x0);
    assert_eq!(setup.mat_1x1.transpose(), setup.mat_1x1);
    assert_eq!(setup.mat_3x3.transpose().transpose(), setup.mat_3x3);

    assert_eq!(Matrix::zeros(2, 3).transpose(), Matrix::zeros(3, 2));
}

#[rstest]
fn submatrix(setup: Fixture) {
    assert_eq!(setup.mat_3x3.submatrix(0, 0), Matrix::from([[5, 6], [8, 9]]));
    assert_eq!(setup.mat_3x3.submatrix(1, 1), Matrix::from([[1, 3], [7, 9]]));
    assert_eq!(setup.mat_3x3.submatrix(0, 2), Matrix::from([[4, 5], [7, 8]]));
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
fn det(setup: Fixture) {
    assert_eq!(setup.mat_0x0.det(), 1.into());
    assert_eq!(setup.mat_1x1.det(), 2.into());
    assert_eq!(setup.mat_3x3.det(), 0.into());

    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).det(), 27.into());
}

#[rstest]
fn inv(setup: Fixture) {
    assert_eq!(setup.mat_0x0.inv(), Some(Matrix::new()));
    assert_eq!(setup.mat_1x1.inv(), Some(Matrix::from([[Fraction::from((1, 2))]])));
    assert_eq!(setup.mat_3x3.inv(), None);

    assert_eq!(
        Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).inv(),
        Some(Matrix::from([
            [Fraction::from((-16, 9)), Fraction::from((8, 9)), Fraction::from((-1, 9))],
            [Fraction::from((14, 9)), Fraction::from((-7, 9)), Fraction::from((2, 9))],
            [Fraction::from((-1, 9)), Fraction::from((2, 9)), Fraction::from((-1, 9))],
        ]))
    );
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
fn lu_decomposition() {
    assert_eq!(
        Matrix::from([[2, 3, 1], [4, 7, 1], [6, 7, 3]]).lu_decomposition(),
        (Matrix::from([[1, 0, 0], [2, 1, 0], [3, -2, 1]]), Matrix::from([[2, 3, 1], [0, 1, -1], [0, 0, -2]]))
    );

    assert_eq!(
        Matrix::from([[2, 3, 1], [4, 0, 1], [6, 7, 3]]).lu_decomposition(),
        (
            Matrix::from([
                [1.into(), 0.into(), 0.into()],
                [2.into(), 1.into(), 0.into()],
                [3.into(), Fraction::from((1, 3)), 1.into()]
            ]),
            Matrix::from([
                [2.into(), 3.into(), 1.into()],
                [0.into(), (-6).into(), (-1).into()],
                [0.into(), 0.into(), Fraction::from((1, 3))]
            ])
        )
    );
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
