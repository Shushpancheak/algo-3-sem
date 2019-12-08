/*
 * C. Пересечение многоугольников
 * Ограничение времени	1 секунда
 * Ограничение памяти	64Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Даны два выпуклых многоугольника на плоскости. В первом n точек, во втором m
 * Определите, пересекаются ли они за O(n + m). Указание. Используйте сумму
 * Минковского.
 */

#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <cmath>

const double EPS = 1e-12;

struct Segment;

struct Vector2 {
  double x;
  double y;

  double GetLength() const;
  double GetLengthSqr() const;
  Vector2& operator+=(const Vector2& vec);
  friend Vector2 operator+(const Vector2& vec_1, const Vector2& vec_2);
  friend Vector2 operator-(const Vector2& vec_1, const Vector2& vec_2);
  friend Vector2 operator-(const Vector2& vec);
  friend Vector2 operator*(double scalar, const Vector2& vec);
  friend double operator*(const Vector2& vec_1, const Vector2& vec_2);
  friend double GetCrossProduct(const Vector2& vec_1, const Vector2& vec_2);
  friend double GetCosOfAngle(const Vector2& vec_1, const Vector2& vec_2); // the least

  ~Vector2() = default;
  Vector2& operator=(const Vector2& vector);
};

struct Point2 {
  Vector2 pos;

  Point2& operator+=(const Vector2& vector);
  friend double GetDistance(const Point2& point_1, const Point2& point_2);
  friend double GetDistanceSqr(const Point2& point_1, const Point2& point_2);
  friend std::istream& operator>>(std::istream& stream, Point2& point);
  friend Point2 operator+(const Point2& point, const Vector2& vec);
  operator Vector2() const {return pos; }

  ~Point2() = default;
};

struct Segment {
  Point2 start;
  Point2 end;

  Vector2 GetDirection() const;
  double GetLength() const;
  double GetLengthSqr() const;
  friend std::istream& operator>>(std::istream& stream, Segment& seg);
  friend Segment operator-(const Point2& point_1, const Point2& point_2);
  friend double GetCosOfAngle(const Segment& seg_1, const Segment& seg_2);
  friend double GetCrossProduct(const Segment& seg_1, const Segment& seg_2);

  ~Segment() = default;
};

struct ConvexPolygon {
  // Clockwise.
  std::vector<Point2> points;

  auto GetLowestRightmostPoint() const;
  auto GetNextPoint(const std::vector<Point2>::const_iterator& point) const;
  friend std::istream& operator>>(std::istream& stream, ConvexPolygon& polygon);
  friend ConvexPolygon operator-(const ConvexPolygon& polygon);
  // Minkowski sums
  friend ConvexPolygon operator+(const ConvexPolygon& polygon_1,
                                 const ConvexPolygon& polygon_2);
  friend ConvexPolygon operator-(const ConvexPolygon& polygon_1,
                                 const ConvexPolygon& polygon_2);
};

bool DoIntersect(const ConvexPolygon& polygon_1, const ConvexPolygon& polygon_2);
bool DoesContainPoint(const ConvexPolygon& polygon, Point2 point);

// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -   - - - MAIN - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -

int main() {
  ConvexPolygon polygon_1, polygon_2;
  std::cin >> polygon_1 >> polygon_2;

  std::cout << (DoIntersect(polygon_1, polygon_2) ? "YES" : "NO") << std::endl;

  return 0;
}


Vector2 operator+(const Vector2& vec_1, const Vector2& vec_2) {
  return {vec_1.x + vec_2.x, vec_1.y + vec_2.y};
}


Vector2 operator-(const Vector2& vec_1, const Vector2& vec_2) {
  return {vec_1.x - vec_2.x, vec_1.y - vec_2.y};
}


Vector2 operator-(const Vector2& vec) {
  return {-vec.x, -vec.y};
}


Vector2 operator*(double scalar, const Vector2& vec) {
  return {scalar * vec.x, scalar * vec.y};
}


double operator*(const Vector2& vec_1, const Vector2& vec_2) {
  return vec_1.x * vec_2.x + vec_1.y * vec_2.y;
}


double GetCrossProduct(const Vector2& vec_1, const Vector2& vec_2) {
  return vec_1.x * vec_2.y - vec_1.y * vec_2.x;
}


double GetCosOfAngle(const Vector2& vec_1, const Vector2& vec_2) {
  return vec_1 * vec_2 / (vec_1.GetLength() * vec_2.GetLength());
}


Segment operator-(const Point2& point_1, const Point2& point_2) {
  Segment seg;
  seg.start = point_2;
  seg.end = point_1;
  return seg;
}


double GetCosOfAngle(const Segment& seg_1, const Segment& seg_2) {
  return GetCosOfAngle(seg_1.GetDirection(), seg_2.GetDirection());
}


double GetCrossProduct(const Segment& seg_1, const Segment& seg_2) {
  return GetCrossProduct(seg_1.GetDirection(), seg_2.GetDirection());
}


std::istream& operator>>(std::istream& stream, ConvexPolygon& polygon) {
  size_t vertices_count = 0;
  std::cin >> vertices_count;
  polygon.points.resize(vertices_count);

  for (size_t i = 0; i < vertices_count; ++i) {
    std::cin >> polygon.points[i];
  }

  return stream;
}


ConvexPolygon operator-(const ConvexPolygon& polygon) {
  ConvexPolygon new_polygon = polygon;

  // Clockwise traversal is preserved.
  for (auto& vertex : new_polygon.points) {
    vertex.pos = -vertex.pos;
  }

  return new_polygon;
}


auto ConvexPolygon::GetLowestRightmostPoint() const {
  auto result_point = points.begin();

  for (auto cur_point = points.begin() + 1; cur_point != points.end(); ++cur_point) {
    if (cur_point->pos.y < result_point->pos.y ||
        cur_point->pos.x > result_point->pos.x &&
        std::abs(cur_point->pos.y - result_point->pos.y) < EPS) {
      result_point = cur_point;
    }
  }

  return result_point;
}


auto ConvexPolygon::GetNextPoint(const std::vector<Point2>::const_iterator& point) const {
  return point + 1 == points.end() ? points.begin() : point + 1;
}


ConvexPolygon operator+(const ConvexPolygon& polygon_1,
    const ConvexPolygon& polygon_2) {
  // Using standard Minkowski sum finding algorithm. 
  const auto start_point_1 = polygon_1.GetLowestRightmostPoint();
  const auto start_point_2 = polygon_2.GetLowestRightmostPoint();
  auto       cur_point_1   = start_point_1;
  auto       cur_point_2   = start_point_2;
  auto       next_point_1  = polygon_1.GetNextPoint(cur_point_1);
  auto       next_point_2  = polygon_2.GetNextPoint(cur_point_2);
  Vector2    cur_dir_1     = (*next_point_1 - *cur_point_1).GetDirection();
  Vector2    cur_dir_2     = (*next_point_2 - *cur_point_2).GetDirection();
  double     cos_1         = 0;
  double     cos_2         = 0;
  Vector2    cur_vec       = {-1, 0}; // The vector which is used in angle comparison.
  unsigned   cur_polygon_to_move = 0; // 1 or 2.

  ConvexPolygon res_polygon;
  res_polygon.points.reserve(polygon_1.points.size() + polygon_2.points.size());
  Point2 cur_res_point = *start_point_1 + *start_point_2;
  res_polygon.points.push_back(cur_res_point);

  do {
    // cur_point == next_point means that all points have been traversed.
    if (cur_point_1 == next_point_1) {
      cur_polygon_to_move = 2;
    } else if (cur_point_2 == next_point_2) {
      cur_polygon_to_move = 1;
    } else {
      cos_1 = GetCosOfAngle(cur_vec, cur_dir_1);
      cos_2 = GetCosOfAngle(cur_vec, cur_dir_2);

      // cos_1 > cos_2 <=> angle_1 < angle_2.
      if (cos_1 > cos_2) {
        cur_polygon_to_move = 1;
      } else {
        cur_polygon_to_move = 2;
      }
    }

    if (cur_polygon_to_move == 1) {
      cur_vec = cur_dir_1;

      if (next_point_1 == start_point_1) {
        cur_point_1 = next_point_1; // leaving next_point the same to imply
                                    // that this is final point for this polygon.
      } else {
        cur_point_1  = next_point_1;
        next_point_1 = polygon_1.GetNextPoint(next_point_1);
        cur_dir_1    = (*next_point_1 - *cur_point_1).GetDirection();
      }
    } else {
      cur_vec = cur_dir_2;

      if (next_point_2 == start_point_2) {
        cur_point_2 = next_point_2; // leaving next_point the same to imply
                                    // that this is final point for this polygon.
      } else {
        cur_point_2  = next_point_2;
        next_point_2 = polygon_2.GetNextPoint(next_point_2);
        cur_dir_2    = (*next_point_2 - *cur_point_2).GetDirection();
      }
    }

    cur_res_point += cur_vec;
    res_polygon.points.push_back(cur_res_point);
  } while (cur_point_1 != start_point_1 || cur_point_2 != start_point_2);

  return res_polygon;
}


ConvexPolygon operator-(const ConvexPolygon& polygon_1,
    const ConvexPolygon& polygon_2) {
  return polygon_1 + (-polygon_2);
}


double GetDistance(const Point2& point_1, const Point2& point_2) {
  return std::sqrt(GetDistanceSqr(point_1, point_2));
}


double GetDistanceSqr(const Point2& point_1, const Point2& point_2) {
  return (point_1.pos.x - point_2.pos.x) * (point_1.pos.x - point_2.pos.x) +
         (point_1.pos.y - point_2.pos.y) * (point_1.pos.y - point_2.pos.y);
}


std::istream& operator>>(std::istream& stream, Point2& point) {
  return stream >> point.pos.x >> point.pos.y;
}


Point2 operator+(const Point2& point, const Vector2& vec) {
  Point2 new_point = point;
  new_point.pos += vec;
  return new_point;
}


Vector2 Segment::GetDirection() const {
  return {end.pos.x - start.pos.x,
          end.pos.y - start.pos.y};
}


double Segment::GetLength() const {
  return GetDistance(start, end);
}


double Segment::GetLengthSqr() const {
  return (start.pos.x - end.pos.x) * (start.pos.x - end.pos.x) +
         (start.pos.y - end.pos.y) * (start.pos.y - end.pos.y);
}


bool DoIntersect(const ConvexPolygon& polygon_1,
    const ConvexPolygon& polygon_2) {
  Point2 zero = {};
  zero.pos = {0, 0};
  return DoesContainPoint(polygon_1 - polygon_2, zero);
}


bool DoesContainPoint(const ConvexPolygon& polygon, Point2 point) {
  // Idea: point should be on the same side for every segment in
  // the polygon.
  // Using cross product to determine to which side the point is on compared to
  // current segment's vector for the sign of the cross product comes from
  // sin of angle between vectors.
  auto start_point = polygon.points.begin();
  auto end_point = polygon.points.begin() + 1;
  int sign = 0;
  int prev_sign = 0;

  for (; end_point != polygon.points.end(); ++end_point, ++start_point) {
    sign =
      std::signbit(GetCrossProduct(*end_point - *start_point, point - *start_point))
        ? -1
        : 1;
    if (prev_sign == 0) {
      prev_sign = sign;
    } else if (prev_sign != sign) {
      return false;
    }
  }

  return true;
}


std::istream& operator>>(std::istream& stream, Segment& seg) {
  return stream >> seg.end >> seg.start;
}


double Vector2::GetLength() const {
  return std::sqrt(GetLengthSqr());
}


double Vector2::GetLengthSqr() const {
  return x * x + y * y;
}


Vector2& Vector2::operator+=(const Vector2& vec) {
  x += vec.x;
  y += vec.y;
  return *this;
}


Vector2& Vector2::operator=(const Vector2& vector) {
  x = vector.x;
  y = vector.y;

  return *this;
}


Point2& Point2::operator+=(const Vector2& vector) {
  pos += vector;
  return *this;
}
