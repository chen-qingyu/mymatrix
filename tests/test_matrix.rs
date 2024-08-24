use mymatrix::Matrix;
use pyinrs::Fraction;
use rstest::{fixture, rstest};

struct Fixture {
    empty: Matrix,
    one: Matrix,
    some: Matrix,
}

#[fixture]
fn setup() -> Fixture {
    Fixture {
        empty: Matrix::new(),
        one: Matrix::from([[1]]),
        some: Matrix::from([[1, 2, 3], [4, 5, 6]]),
    }
}

#[rstest]
fn basics(setup: Fixture) {
    assert_eq!(setup.empty.row_size(), 0);
    assert_eq!(setup.empty.col_size(), 0);
    assert!(setup.empty.is_empty());

    assert_eq!(setup.one.row_size(), 1);
    assert_eq!(setup.one.col_size(), 1);
    assert!(!setup.one.is_empty());

    assert_eq!(setup.some.row_size(), 2);
    assert_eq!(setup.some.col_size(), 3);
    assert!(!setup.some.is_empty());
}

#[rstest]
fn compare(setup: Fixture) {
    assert!(setup.some == Matrix::from([[1, 2, 3], [4, 5, 6]]));
    assert!(setup.some != Matrix::from([[1, 2, 3], [4, 5, 0]]));
}

#[rstest]
fn access(mut setup: Fixture) {
    for r in 0..setup.some.row_size() {
        for c in 0..setup.some.col_size() {
            assert_eq!(setup.some[r][c], (r as i32 * 3 + c as i32 + 1).into());
        }
    }

    setup.some[0][0] = 0.into();
    assert_eq!(setup.some[0][0], 0.into());
}

#[rstest]
fn rank() {
    assert_eq!(Matrix::ones(2, 2).rank(), 1);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).rank(), 3);
}

#[rstest]
fn det() {
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).det(), Some(0.into()));
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).det(), Some(27.into()));
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).det(), None);
}

#[rstest]
fn inv() {
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]).inv(),
        Some(Matrix::from([
            [Fraction::from(-2), Fraction::from(1)],
            [Fraction::from((3, 2)), Fraction::from((-1, 2))]
        ]))
    );

    assert_eq!(
        Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).inv(),
        Some(Matrix::from([
            [Fraction::from((-16, 9)), Fraction::from((8, 9)), Fraction::from((-1, 9))],
            [Fraction::from((14, 9)), Fraction::from((-7, 9)), Fraction::from((2, 9))],
            [Fraction::from((-1, 9)), Fraction::from((2, 9)), Fraction::from((-1, 9))],
        ]))
    );

    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).inv(), None);
}

#[rstest]
fn submatrix() {
    let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    assert_eq!(m.submatrix(0, 0), Matrix::from([[5, 6], [8, 9]]));
    assert_eq!(m.submatrix(1, 1), Matrix::from([[1, 3], [7, 9]]));
    assert_eq!(m.submatrix(0, 2), Matrix::from([[4, 5], [7, 8]]));
}

#[rstest]
fn minor() {
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).minor(), Matrix::from([[4, 3], [2, 1]]));
    assert_eq!(
        Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).minor(),
        Matrix::from([[-3, -6, -3], [-6, -12, -6], [-3, -6, -3]])
    );
}

#[rstest]
fn cofactor() {
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).cofactor(), Matrix::from([[4, -3], [-2, 1]]));
    assert_eq!(
        Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).cofactor(),
        Matrix::from([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])
    );
}

#[rstest]
fn adj() {
    let m1 = Matrix::from([[1, 2], [3, 4]]);
    assert_eq!(m1.adj(), Matrix::from([[4, -2], [-3, 1]]));
    assert_eq!(m1.adj() * (Fraction::from(1) / m1.det().unwrap()), m1.inv().unwrap()); // inv(A) = (1/det(A)) * adj(A)

    let m2 = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]);
    assert_eq!(m2.adj(), Matrix::from([[-48, 24, -3], [42, -21, 6], [-3, 6, -3]]));
    assert_eq!(m2.adj() * (Fraction::from(1) / m2.det().unwrap()), m2.inv().unwrap());
}

#[rstest]
fn is_symmetric(setup: Fixture) {
    assert_eq!(setup.empty.is_symmetric(), true);
    assert_eq!(setup.one.is_symmetric(), true);
    assert_eq!(setup.some.is_symmetric(), false);
    assert_eq!((&setup.some * &setup.some.transpose()).is_symmetric(), true);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).is_symmetric(), false);
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
fn elementary_row_operations() {
    let mut matrix = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    assert_eq!(matrix.e_row_swap(0, 1), &Matrix::from([[4, 5, 6], [1, 2, 3], [7, 8, 9]]));
    assert_eq!(matrix.e_scalar_multiplication(1, 2.into()), &Matrix::from([[4, 5, 6], [2, 4, 6], [7, 8, 9]]));
    assert_eq!(matrix.e_row_sum(0, 1, (-1).into()), &Matrix::from([[2, 1, 0], [2, 4, 6], [7, 8, 9]]));
}

#[rstest]
fn rref() {
    assert_eq!(Matrix::ones(2, 2).rref(), Matrix::from([[1, 1], [0, 0]]));

    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).rref(), Matrix::from([[1, 0, -1], [0, 1, 2]]));
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).rref(), Matrix::from([[1, 0], [0, 1], [0, 0]]));

    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6], [7, 8]]).rref(), Matrix::from([[1, 0], [0, 1], [0, 0], [0, 0]]));
    assert_eq!(Matrix::from([[1, 2, 3, 4], [5, 6, 7, 8]]).rref(), Matrix::from([[1, 0, -1, -2], [0, 1, 2, 3]]));
}

#[rstest]
fn lu_decomposition() {
    assert_eq!(
        Matrix::from([[2, 3, 1], [4, 7, 1], [6, 7, 3]]).lu_decomposition(),
        (Matrix::from([[1, 0, 0], [2, 1, 0], [3, -2, 1]]), Matrix::from([[2, 3, 1], [0, 1, -1], [0, 0, -2]]))
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
fn transpose() {
    assert_eq!(Matrix::create(2, 3, 1.into()).transpose(), Matrix::create(3, 2, 1.into()));
    assert_eq!(Matrix::create(1, 3, 3.into()).transpose(), Matrix::create(3, 1, 3.into()));
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
    assert_eq!(format!("{}", setup.empty), "[\n]");
    assert_eq!(format!("{}", setup.one), "[\n1\n]");
    assert_eq!(format!("{}", setup.some), "[\n1 2 3\n4 5 6\n]");

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
