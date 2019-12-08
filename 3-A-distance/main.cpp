/*
 * A. Расстояние между отрезками
 * Ограничение времени	1 секунда
 * Ограничение памяти	64Mb
 * Ввод	стандартный ввод
 * Вывод	стандартный вывод
 * Даны два отрезка в пространстве (x1, y1, z1) - (x2, y2, z2) и (x3, y3, z3) - (x4,
 * y4, z4). Найдите расстояние между отрезками.
 *
 * Формат ввода
 * Заданы целые x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4. Координаты по модулю
 * не превосходят 1000.
 *
 * Формат вывода
 * Выведите искомую величину c точностью не менее 6 знаков после десятичной точки.
 */

#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <algorithm>
#include <cmath>

const double EPS = 1e-12;

struct Segment;

struct Vector3 {
  double x;
  double y;
  double z;

  double GetLength() const;
  double GetLengthSqr() const;
  Vector3& operator+=(const Vector3& vec);
  friend Vector3 operator+(const Vector3& vec_1, const Vector3& vec_2);
  friend Vector3 operator-(const Vector3& vec_1, const Vector3& vec_2);
  friend Vector3 operator-(const Vector3& vec);
  friend Vector3 operator*(double scalar, const Vector3& vec);
  friend double operator*(const Vector3& vec_1, const Vector3& vec_2);

  ~Vector3() = default;
  Vector3& operator=(const Vector3& vector3);
};

struct Point3 {
  Vector3 pos;

  friend double GetDistance(const Point3& point_1, const Point3& point_2);
  friend double GetDistanceSqr(const Point3& point_1, const Point3& point_2);
  friend std::istream& operator>>(std::istream& stream, Point3& point);
  friend Point3 operator+(const Point3& point, const Vector3 vec);
  operator Vector3() const {return pos; };

  ~Point3() = default;
};

struct Line {
  virtual ~Line() = default;
};

struct LineEquation : Line {
  double A;
  double B;
  double C;

  ~LineEquation() = default;
};

struct LineDirection : Line {
  Point3 point;
  Vector3 direction;

  friend Segment GetDistanceSegment(const LineDirection& line_1,
                                    const LineDirection& line_2);
  friend std::pair<double, double> GetDistanceSegmentScalars(
      const LineDirection& line_1,
      const LineDirection& line_2);

  ~LineDirection() = default;
};

struct Segment {
  Point3 start;
  Point3 end;

  Vector3 GetDirection() const;
  double GetLength() const;
  double GetLengthSqr() const;
  friend std::istream& operator>>(std::istream& stream, Segment& seg);
  friend double GetDistance(const Segment& seg_1, const Segment& seg_2);
  friend double GetDistanceSqr(const Segment& seg_1, const Segment& seg_2);
  friend Segment operator-(const Point3& point_1, const Point3& point_2);
  explicit operator LineDirection() const;

  ~Segment() = default;
};

Segment GetParameterizedBridgeBetweenLines(
    const LineDirection& line_1,
    const LineDirection& line_2,
    double scalar_1,
    double scalar_2);

// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -   - - - MAIN - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -

int main() {
  Segment seg_1 = {}, seg_2 = {};

  std::cin >> seg_1 >> seg_2;

  std::cout.precision(12);
  std::cout << std::fixed << GetDistance(seg_1, seg_2);

  return 0;
}


Vector3 operator+(const Vector3& vec_1, const Vector3& vec_2) {
  return {vec_1.x + vec_2.x, vec_1.y + vec_2.y, vec_1.z + vec_2.z};
}


Vector3 operator-(const Vector3& vec_1, const Vector3& vec_2) {
  return {vec_1.x - vec_2.x, vec_1.y - vec_2.y, vec_1.z - vec_2.z};
}


Vector3 operator-(const Vector3& vec) {
  return {-vec.x, -vec.y, -vec.z};
}


Vector3 operator*(double scalar, const Vector3& vec) {
  return {scalar * vec.x, scalar * vec.y, scalar * vec.z};
}


double operator*(const Vector3& vec_1, const Vector3& vec_2) {
  return vec_1.x * vec_2.x + vec_1.y * vec_2.y + vec_1.z * vec_2.z;
}


Segment operator-(const Point3& point_1, const Point3& point_2) {
  Segment seg;
  seg.start = point_2;
  seg.end = point_1;
  return seg;
}


double GetDistance(const Point3& point_1, const Point3& point_2) {
  return std::sqrt(GetDistanceSqr(point_1, point_2));
}


double GetDistanceSqr(const Point3& point_1, const Point3& point_2) {
  return (point_1.pos.x - point_2.pos.x) * (point_1.pos.x - point_2.pos.x) +
         (point_1.pos.y - point_2.pos.y) * (point_1.pos.y - point_2.pos.y) +
         (point_1.pos.z - point_2.pos.z) * (point_1.pos.z - point_2.pos.z);
}


std::istream& operator>>(std::istream& stream, Point3& point) {
  return stream >> point.pos.x >> point.pos.y >> point.pos.z;
}


Point3 operator+(const Point3& point, const Vector3 vec) {
  Point3 new_point = point;
  new_point.pos += vec;
  return new_point;
}


Vector3 Segment::GetDirection() const {
  return {end.pos.x - start.pos.x,
          end.pos.y - start.pos.y,
          end.pos.z - start.pos.z};
}


double Segment::GetLength() const {
  return GetDistance(start, end);
}


double Segment::GetLengthSqr() const {
  return (start.pos.x - end.pos.x) * (start.pos.x - end.pos.x) +
         (start.pos.y - end.pos.y) * (start.pos.y - end.pos.y) +
         (start.pos.z - end.pos.z) * (start.pos.z - end.pos.z);
}


Segment::operator LineDirection() const {
  LineDirection line;
  line.point = start;
  line.direction = (end - start).GetDirection();
  return line;
}


Segment GetParameterizedBridgeBetweenLines(
    const LineDirection& line_1,
    const LineDirection& line_2,
    double scalar_1,
    double scalar_2) {
  Segment seg;
  seg.start = line_1.point + scalar_1 * line_1.direction;
  seg.end = line_2.point + scalar_2 * line_2.direction;
  return seg;
}


Segment GetDistanceSegment(
    const LineDirection& line_1,
    const LineDirection& line_2) {
  // See GetDistanceSegmentScalars for explanation.
  auto [s_c, t_c] = GetDistanceSegmentScalars(line_1, line_2);
  return GetParameterizedBridgeBetweenLines(line_1, line_2, s_c, t_c);
}


std::pair<double, double> GetDistanceSegmentScalars(const LineDirection& line_1,
    const LineDirection& line_2) {
  // Explanation:
  // (See segments distance first)
  //
  // Trying to find (s_c, t_c) in vec(s, t) = seg_1(s) - seg_2(t) =
  // = vec_0 + s * dir_1 - t * dir_2.
  //
  // We know that vec(s_c, t_c) is perpendicular to dir_1 and dir_2,
  // which is enough for a solvable equation:
  // { (dir_1 * dir_1)s_c - (dir_1 * dir_2)t_c = dir_1 * vec_0,
  // { (dir_2 * dir_1)s_c - (dir_2 * dir_2)t_c = dir_2 * vec_0.
  //
  // Let a = dir_1 * dir_1, b = dir_1 * dir_2, c = dir_2 * dir_2,
  // d = dir_1 * vec_0, e = dir_2 * vec_0;
  // => (s_c, t_c) = ( (be - cd) / (ac - bb), (ae - bd) / (ac - bb) ).
  //
  // If ac - bb is nearly zero, suppose it's zero (which means that lines are
  // parallel) and fixate s_c = 0, then we get:
  // t_c = d / b = e / c;
  //
  // Note that ac - bb = |dir_1|^2|dir_2|^2 - |dir_1|^2|dir_2|^2 cos^2(phi) =
  // = (|dir_1||dir_2|sin(phi))^2 = cross_product.

  const Vector3 dir_1 = line_1.direction;
  const Vector3 dir_2 = line_2.direction;
  const Vector3 vec_0 = (line_1.point - line_2.point).GetDirection();

  const double a = dir_1 * dir_1;
  const double b = dir_1 * dir_2;
  const double c = dir_2 * dir_2;
  const double d = dir_1 * vec_0;
  const double e = dir_2 * vec_0;
  const double cross_product = a * c - b * b;

  double s_c = 0, t_c = 0;

  // If parallel or one of these is dot.
  if (std::abs(cross_product) < EPS) { 
    s_c = 0;
    t_c = (b > c ? d / b : e / c); // Use the largest denominator for preciseness.
  } else {
    s_c = (b * e - c * d) / cross_product;
    t_c = (a * e - b * d) / cross_product;
  }

  return std::make_pair(s_c, t_c);
}


std::istream& operator>>(std::istream& stream, Segment& seg) {
  return stream >> seg.end >> seg.start;
}


double GetDistance(const Segment& seg_1, const Segment& seg_2) {
  return std::sqrt(GetDistanceSqr(seg_1, seg_2));
}


double GetDistanceSqr(const Segment& seg_1, const Segment& seg_2) {
  // Explanation:
  // seg_i    = [seg_i.start, seg_i.end];
  // seg_i(t) = seg_i.start + (seg_i.end - seg_i.start) * t = seg_i.start + dir_i * t;
  // (t \in [0; 1]).
  //
  // Let vec(s, t) = seg_1(s) - seg_2(t) => We need to minimize vec(s, t);
  // (s, t \in [0; 1]).
  //
  // => Minimize |vec|^2 = vec * vec =
  // = (vec_0 + s * dir_1 - t * dir_2) * (vec_0 + s * dir_1 - t * dir_2),
  // where vec_0 = seg_1.start - seg_2.start;
  //
  // vec^2 defines a paraboloid in (s, t) plane with minimum at C = (s_c, t_c);
  // We are restricted to a region G = {(s, t) | s, t \in [0; 1]}.
  //
  // For Lines:
  // Trying to find (s_c, t_c) in vec(s, t) = seg_1(s) - seg_2(t) =
  // = vec_0 + s * dir_1 - t * dir_2.
  //
  // We know that vec(s_c, t_c) is perpendicular to dir_1 and dir_2,
  // which is enough for a solvable equation:
  // { (dir_1 * dir_1)s_c - (dir_1 * dir_2)t_c = dir_1 * vec_0,
  // { (dir_2 * dir_1)s_c - (dir_2 * dir_2)t_c = dir_2 * vec_0.
  //
  // Let a = dir_1 * dir_1, b = dir_1 * dir_2, c = dir_2 * dir_2,
  // d = dir_1 * vec_0, e = dir_2 * vec_0;
  // => (s_c, t_c) = ( (be - cd) / (ac - bb), (ae - bd) / (ac - bb) ),
  // if resulting segment has both of his ends on given segments.
  //
  // Algorithm is as follows:
  // 1. Find min vector between Lines on which segments reside.
  //    If (s_c, s_t) is in G, we found distance.
  //
  // 2. If not, We will search on edges of square G:
  //    if we are on either edge, then one of s, t is 0 or 1. That leaves only
  //    one argument in function vec and makes it easily differentiable.
  //    Suppose s == 0, then we can find argmin t_0 from:
  //       0   = d/dt(|vec|^2) = -2 * dir_2 * (vec_0 - t_0 * dir_2) =>
  //    => t_0 = (vec_0 * dir_2) / (dir_2 * dir_2).
  //
  //    If t_0 is within [0; 1], then we found points, otherwise --
  //    search on other sides.

  const LineDirection seg_line_1 = static_cast<LineDirection>(seg_1);
  const LineDirection seg_line_2 = static_cast<LineDirection>(seg_2);

  auto [s_c, t_c] =
      GetDistanceSegmentScalars(seg_line_1, seg_line_2);

  if (0 <= s_c && s_c <= 1 &&
      0 <= t_c && t_c <= 1) {
    return GetParameterizedBridgeBetweenLines(
              seg_line_1, seg_line_2, s_c, t_c).GetLengthSqr();
  }

  const Vector3 dir_1 = seg_1.GetDirection();
  const Vector3 dir_2 = seg_2.GetDirection();
  const Vector3 vec_0 = (seg_1.start - seg_2.start).GetDirection();
  const double a = dir_1 * dir_1; // always >= 0
  const double b = dir_1 * dir_2;
  const double c = dir_2 * dir_2; // always >= 0
  const double d = dir_1 * vec_0;
  const double e = dir_2 * vec_0;
  const double cross_product = a * c - b * b; // always >= 0

  // Non parallel segments.
  if (cross_product > 0) {
    // (s_c, t_c) = ( (be - cd) / (ac - bb), (ae - bd) / (ac - bb) ).
    const double be = b * e;
    const double cd = c * d;

    // s_c <= 0
    if (be <= cd) {
      // e <= 0 -> ae <= 0 -> t_c <= 0
      if (e <= 0) {
        s_c = (-d < a ? (-d > 0 ? -d / a : 0) : 1);
        t_c = 0;

      // 0 < t_c < 1
      } else if (e < c) {
        s_c = 0;
        t_c = e / c;

      // t_c >= 1
      } else {
        s_c = (b - d < a ? (b - d > 0 ? (b - d) / a : 0) : 1);
        t_c = 1;
      }

    // s_c > 0
    } else {
      s_c = be - cd;

      // s_c >= 1
      if (s_c >= cross_product) {
        // t_c <= 0
        if (b + e <= 0) {
          s_c = (-d > 0 ? (-d < a ? -d / a : 1) : 0);
          t_c = 0;

        // 0 < t_c < 1
        } else if (b + e < c) {
          s_c = 1;
          t_c = (b + e) / c;

        // t_c >= 1
        } else {
          s_c = b - d > 0 ? (b - d < a ? (b - d) / a : 1) : 0;
          t_c = 1;
        }

      // 0 < s_c < 1
      } else {
        const double ae = a * e;
        const double bd = b * d;

        // t_c <= 0
        if (ae <= bd) {
          s_c = -d > 0 ? (-d < a ? -d / a : 1) : 0;
          t_c = 0;

        // t_c > 0
        } else {
          t_c = ae - bd;

          // t_c >= 1
          if (t_c >= cross_product) {
            s_c = b - d > 0 ? (b - d < a ? (b - d) / a : 1) : 0;
            t_c = 1;

          // 0 < t_c < 1
          } else {
            s_c /= cross_product;
            t_c /= cross_product;
          }
        }
      }
    }

  // Parallel segments
  } else {

    // t_c <= 0 (setting to EPS to also catch conditions when segments are points).
    if (e <= EPS) {
      s_c = -d > 0 ? (-d < a ? -d / a : 1) : 0;
      t_c = 0;

    // t_c >= 1
    } else if (e >= c) {
      s_c = (b - d > 0 ? (b - d < a ? (b - d) / a : 1) : 0);
      t_c = 1;

    // s_c <= 0
    } else {
      s_c = 0;
      t_c = e / c;
    }
  }

  return GetParameterizedBridgeBetweenLines(
              seg_line_1, seg_line_2, s_c, t_c).GetLengthSqr();
}


double Vector3::GetLength() const {
  return std::sqrt(GetLengthSqr());
}


double Vector3::GetLengthSqr() const {
  return x * x + y * y + z * z;
}


Vector3& Vector3::operator+=(const Vector3& vec) {
  x += vec.x;
  y += vec.y;
  z += vec.z;
  return *this;
}


Vector3& Vector3::operator=(const Vector3& vector3) {
  x = vector3.x;
  y = vector3.y;
  z = vector3.z;

  return *this;
}
