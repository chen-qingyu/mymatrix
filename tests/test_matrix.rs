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
    assert_eq!(Matrix::create(2, 2, 1.into()).rank(), 1);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rank(), 2);
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).rank(), 3);
}

#[rstest]
fn det() {
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).det(), 0.into());
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).det(), 27.into());
}

#[rstest]
#[should_panic(expected = "Error: The dimensions mismatch.")]
fn bad_det() {
    Matrix::from([[1, 2, 3], [4, 5, 6]]).det();
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
}

#[rstest]
fn append() {
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]).append_row(Matrix::create(2, 2, 0.into())),
        &Matrix::from([[1, 2], [3, 4], [0, 0], [0, 0]])
    );
    assert_eq!(
        Matrix::from([[1, 2, 3, 4, 5]]).append_row(Matrix::from([[6, 7, 8, 9, 10]])),
        &Matrix::from([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
    );

    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]).append_col(Matrix::create(2, 2, 0.into())),
        &Matrix::from([[1, 2, 0, 0], [3, 4, 0, 0]])
    );
    assert_eq!(
        Matrix::from([[1, 2, 3, 4, 5]]).append_col(Matrix::from([[6, 7, 8, 9, 10]])),
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
fn to_row_echelon() {
    assert_eq!(Matrix::create(2, 2, 1.into()).to_row_echelon(), &Matrix::from([[1, 1], [0, 0]]));
    assert_eq!(Matrix::from([[1, 2, 3], [4, 5, 6]]).to_row_echelon(), &Matrix::from([[1, 2, 3], [0, -3, -6]]));
    assert_eq!(Matrix::from([[1, 2], [3, 4], [5, 6]]).to_row_echelon(), &Matrix::from([[1, 2], [0, -2], [0, 0]]));
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
    assert_eq!(Matrix::create(2, 3, 1.into()) * Matrix::create(2, 3, 2.into()), Matrix::create(2, 3, 2.into()));
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]) * Matrix::from([[-1, -1], [-1, -1]]),
        Matrix::from([[-1, -2], [-3, -4]])
    );

    assert_eq!(Matrix::create(2, 3, 1.into()) * Fraction::from(2), Matrix::create(2, 3, 2.into()));
    assert_eq!(Matrix::from([[1, 2], [3, 4]]) * Fraction::from(3), Matrix::from([[3, 6], [9, 12]]));
}

#[rstest]
fn dot() {
    assert_eq!(
        Matrix::dot(&Matrix::create(2, 2, 1.into()), &Matrix::create(2, 2, 2.into())),
        Matrix::create(2, 2, 4.into())
    );
    assert_eq!(
        Matrix::dot(&Matrix::create(1, 3, 1.into()), &Matrix::create(3, 1, 1.into())),
        Matrix::create(1, 1, 3.into())
    );
}

#[rstest]
fn format(setup: Fixture) {
    assert_eq!(format!("{}", setup.empty), "[[]]");
    assert_eq!(format!("{}", setup.one), "[[\n1;\n]]");
    assert_eq!(format!("{}", setup.some), "[[\n1 2 3;\n4 5 6;\n]]");
}
