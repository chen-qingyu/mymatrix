use mymatrix::Vector;
use pyinrs::Fraction;
use rstest::{fixture, rstest};

struct Fixture {
    empty: Vector,
    one: Vector,
    some: Vector,
}

#[fixture]
fn setup() -> Fixture {
    Fixture {
        empty: Vector::new(),
        one: Vector::from([1]),
        some: Vector::from([1, 2, 3, 4, 5]),
    }
}

#[rstest]
fn basics(setup: Fixture) {
    assert_eq!(setup.empty.size(), 0);
    assert!(setup.empty.is_empty());

    assert_eq!(setup.one.size(), 1);
    assert!(!setup.one.is_empty());

    assert_eq!(setup.some.size(), 5);
    assert!(!setup.some.is_empty());
}

#[rstest]
fn compare(setup: Fixture) {
    assert!(setup.some == Vector::from([1, 2, 3, 4, 5]));
    assert!(setup.some != Vector::from([1, 3, 5]));
}

#[rstest]
fn access(mut setup: Fixture) {
    // for
    for i in 0..setup.some.size() {
        assert_eq!(setup.some[i], ((i + 1) as i32).into());
    }

    // assignment
    setup.some[0] = 0.into();
    assert_eq!(setup.some[0], 0.into());
}

#[rstest]
fn is_zero() {
    assert!(Vector::from([0]).is_zero());
    assert!(Vector::from([0, 0, 0]).is_zero());

    assert!(!Vector::from([1]).is_zero());
    assert!(!Vector::from([0, 0, 1]).is_zero());
}

#[rstest]
fn is_orthogonal() {
    let zero = Vector::from([0, 0]);
    assert!(zero.is_orthogonal(&Vector::from([0, 0])));
    assert!(zero.is_orthogonal(&Vector::from([1, 1])));
    assert!(zero.is_orthogonal(&Vector::from([2, 3])));

    let one = Vector::from([1, 1]);
    assert!(!one.is_orthogonal(&Vector::from([1, 1])));
    assert!(one.is_orthogonal(&Vector::from([1, -1])));
    assert!(one.is_orthogonal(&Vector::from([-1, 1])));
    assert!(one.is_orthogonal(&Vector::from([-2, 2])));
}

#[rstest]
fn is_parallel() {
    let zero = Vector::from([0, 0]);
    assert!(zero.is_parallel(&Vector::from([0, 0])));
    assert!(zero.is_parallel(&Vector::from([1, 1])));
    assert!(zero.is_parallel(&Vector::from([2, 3])));

    let some = Vector::from([3, 4]);
    assert!(!some.is_parallel(&Vector::from([1, 1])));
    assert!(some.is_parallel(&Vector::from([3, 4])));
    assert!(some.is_parallel(&Vector::from([-3, -4])));
    assert!(some.is_parallel(&Vector::from([6, 8])));

    assert!(Vector::from([1, 1]).is_parallel(&Vector::from([2, 2])));
}

#[rstest]
fn norm() {
    assert_eq!(Vector::from([0]).norm(), 0.0);
    assert_eq!(Vector::from([1]).norm(), 1.0);
    assert_eq!(Vector::from([3, 4]).norm(), 5.0);
}

#[rstest]
fn count_leading_zeros() {
    assert_eq!(Vector::from([0]).count_leading_zeros(), 1);
    assert_eq!(Vector::from([0, 1]).count_leading_zeros(), 1);
    assert_eq!(Vector::from([0, 0]).count_leading_zeros(), 2);
    assert_eq!(Vector::from([0, 0, 1]).count_leading_zeros(), 2);
    assert_eq!(Vector::from([0, 0, 0, 1, 2, 3]).count_leading_zeros(), 3);
    assert_eq!(Vector::from([0, 0, 0, 1, 2, 3, 0, 0, 0]).count_leading_zeros(), 3);
}

#[rstest]
fn cross() {
    assert_eq!(Vector::cross(&[1, 2].into(), &[3, 4].into()), [-2].into());
    assert_eq!(Vector::cross(&[1, 2, 3].into(), &[4, 5, 6].into()), [-3, 6, -3].into());
}

#[rstest]
#[should_panic(expected = "Error: Incompatible dimensions for cross product.")]
fn bad_cross() {
    Vector::cross(&[1].into(), &[2].into());
}

#[rstest]
fn add() {
    assert_eq!(Vector::from([1]) + Vector::from([1]), Vector::from([2]));
    assert_eq!(Vector::from([1, 2, 3]) + Vector::from([4, 5, 6]), Vector::from([5, 7, 9]));
    assert_eq!(Vector::from([1, 2]) + Vector::from([1, 2]) + Vector::from([1, 2]), Vector::from([3, 6]));
}

#[rstest]
#[should_panic(expected = "Error: The dimensions mismatch.")]
fn bad_add() {
    let _ = Vector::from([1]) + Vector::from([1, 2]);
}

#[rstest]
fn sub() {
    assert_eq!(Vector::from([1]) - Vector::from([1]), Vector::from([0]));
    assert_eq!(Vector::from([1, 2, 3]) - Vector::from([4, 5, 6]), Vector::from([-3, -3, -3]));
    assert_eq!(Vector::from([1, 2]) - Vector::from([1, 2]) - Vector::from([1, 2]), Vector::from([-1, -2]));
}

#[rstest]
fn scalar_mul() {
    assert_eq!(Fraction::from(1) * Vector::from([1]), Vector::from([1]));
    assert_eq!(Vector::from([1]) * Fraction::from(1), Vector::from([1]));
    assert_eq!(Vector::from([1, 2, 3]) * Fraction::from(2), Vector::from([2, 4, 6]));
    assert_eq!(
        Vector::from([1, 2]) * Fraction::from(2) * Fraction::from((4, 10)),
        Vector::from([Fraction::from((8, 10)), Fraction::from((16, 10))])
    );
}

#[rstest]
fn dot_mul() {
    assert_eq!(Vector::from([1]) * Vector::from([1]), 1.into());
    assert_eq!(Vector::from([1, 2, 3]) * Vector::from([4, 5, 6]), 32.into());
}

#[rstest]
fn format(setup: Fixture) {
    assert_eq!(format!("{}", setup.empty), "[]");
    assert_eq!(format!("{}", setup.one), "[1]");
    assert_eq!(format!("{}", setup.some), "[1 2 3 4 5]");

    assert_eq!(
        format!("{}", Vector::from([Fraction::from((-3, 4)), Fraction::new(), Fraction::from((5, 6))])),
        "[-3/4    0  5/6]"
    );
}
