/*
 * B. Выпуклая оболочка 3D
 * Ограничение времени	0.2 секунды
 * Ограничение памяти	1Mb
 * Ввод	стандартный ввод или input.txt
 * Вывод	стандартный вывод или output.txt
 * Даны n точек в пространстве. Никакие 4 точки не лежат в одной плоскости. Найдите
 * выпуклую оболочку этих точек за O(n*log(n)).
 *
 * Формат ввода
 * Первая строка содержит число m - количество тестов. В последующих строках описаны
 * сами тесты. Каждый тест начинается со строки, содержащей n (n ≤ 1000) - число
 * точек. Далее, в n строках даны по три числа - координаты точек. Все координаты
 * целые, не превосходят по модулю 500.
 *
 * Формат вывода
 * Для каждого теста выведите следующее. В первую строку выведите количество граней m.
 * Далее в последующие m строк выведите описание граней: количество точек грани (=3) и
 * номера точек в исходном множестве. Точки нумеруются в том же порядке, в котором они
 * даны во входном файле. Точки в пределах грани должны быть отсортированы в порядке
 * против часовой стрелки относительно внешней нормали к грани. Первая точка – точка с
 * минимальным номером. Порядок граней лексикографический.
 *
 * Chan's algorithm.
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _MSC_VER
#include <basetsd.h>
typedef SSIZE_T ssize_t;
#endif

const double INF          = 1e99;
const double EPS          = 1e-12;
const double EPS_ROTATION = 1e-2;

struct Segment;

struct Vector3 {
  double x = 0;
  double y = 0;
  double z = 0;

  double GetLength() const;
  double GetLengthSqr() const;
  void   RotateOX(double radians);
  void   RotateOY(double radians);
  void   RotateOZ(double radians);

  Vector3&       operator+=(const Vector3& vec);
  friend Vector3 operator+(const Vector3& vec_1, const Vector3& vec_2);
  friend Vector3 operator-(const Vector3& vec_1, const Vector3& vec_2);
  friend Vector3 operator-(const Vector3& vec);
  friend Vector3 operator*(double scalar, const Vector3& vec);
  friend double  operator*(const Vector3& vec_1, const Vector3& vec_2);
  friend double  GetCrossProductXY(const Vector3& vec_1, const Vector3& vec_2);
  friend double  GetCrossProductXZ(const Vector3& vec_1, const Vector3& vec_2);

  ~Vector3() = default;
  Vector3& operator=(const Vector3& vector3);
};

struct Point3 {
  Vector3 pos = {};

  friend double GetDistance(const Point3& point_1, const Point3& point_2);
  friend double GetDistanceSqr(const Point3& point_1, const Point3& point_2);
  friend std::istream& operator>>(std::istream& stream, Point3& point);
  friend Point3 operator+(const Point3& point, const Vector3 vec);
  operator Vector3() const { return pos; }

  ~Point3() = default;
};

struct LinkedPoint3 : Point3 {
  size_t id = 0;
  LinkedPoint3* next_point = nullptr;
  LinkedPoint3* prev_point = nullptr;

  // @return true if it was insertion, false - deletion.
  bool DeleteOrInsert();
  friend bool operator<(const LinkedPoint3& p1, const LinkedPoint3& p2);
};

typedef LinkedPoint3 PointsEvent;

struct Segment {
  Point3 start = {};
  Point3 end   = {};

  Vector3 GetDirection() const;
  double GetLength() const;
  double GetLengthSqr() const;
  friend Segment operator-(const Point3& point_1, const Point3& point_2);

  ~Segment() = default;
};

struct Face {
  Face(std::initializer_list<size_t> list);

  size_t point_id[3];

  void MakeMinIdPointFirst();
  friend bool operator<(const Face& f1, const Face& f2);
};

struct Hull {
  std::vector<Face> faces;

  void SortFacesLexicographically();
  friend std::ostream& operator<<(std::ostream& stream, const Hull& hull);
};


double GetCrossProductOfPoints(
    const LinkedPoint3* point_1,
    const LinkedPoint3* point_2,
    const LinkedPoint3* point_3);

double GetTimeOfTurning(
    const LinkedPoint3* point_1,
    const LinkedPoint3* point_2,
    const LinkedPoint3* point_3);

std::vector<PointsEvent*> GetEventsInHullBuilding(
    std::vector<LinkedPoint3>& points,
    size_t begin,
    size_t end);

Hull BuildConvexHull(std::vector<LinkedPoint3> points);

// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -   - - - MAIN - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -

int main() {
  size_t tests_count  = 0;
  size_t points_count = 0;

  std::cin >> tests_count;

  for (size_t i = 0; i < tests_count; ++i) {
    std::cin >> points_count;
    std::vector<LinkedPoint3> points;
    points.reserve(points_count);

    for (size_t j = 0; j < points_count; j++) {
      LinkedPoint3 point;
      std::cin >> point;
      point.id = j;
      point.pos.RotateOX(EPS_ROTATION);
      point.pos.RotateOY(EPS_ROTATION);
      point.pos.RotateOZ(EPS_ROTATION);
      points.push_back(point);
    }

    Hull convex_hull = BuildConvexHull(points);
    convex_hull.SortFacesLexicographically();

    std::cout << convex_hull;
  }
  return 0;
}


// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - - IMPLEMENTATIONS - - - - - - - - - - -  - - - - - - - -
// - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - - -  - - - - - - - -

std::pair<LinkedPoint3*, LinkedPoint3*> FindBridge(
    std::vector<LinkedPoint3>& sorted_points,
    const size_t mid) {

  LinkedPoint3* bridge_point_1 = &sorted_points[mid - 1];
  LinkedPoint3* bridge_point_2 = &sorted_points[mid];

  // Moving points until both u- u v and u v v+ make counter-clockwise rotation.
  for (; ; ) {
    if (GetCrossProductOfPoints(bridge_point_1, bridge_point_2,
                                bridge_point_2->next_point) < 0) {
      bridge_point_2 = bridge_point_2->next_point;
    } else if (GetCrossProductOfPoints(bridge_point_1->prev_point,
                                       bridge_point_1, bridge_point_2) < 0) {
      bridge_point_1 = bridge_point_1->prev_point;
    } else {
      break;
    }
  }

  return std::make_pair(bridge_point_1, bridge_point_2);
}


void ReverseDeletionAndInsertion(
    std::vector<LinkedPoint3>& points,
    const size_t middle,
    std::vector<PointsEvent*>& res_events,
    LinkedPoint3* bridge_pt_1,
    LinkedPoint3* bridge_pt_2) {
  // Going back in time to update pointers.
  bridge_pt_1->next_point = bridge_pt_2;
  bridge_pt_2->prev_point = bridge_pt_1;
  for (ssize_t i = static_cast<ssize_t>(res_events.size() - 1); i >= 0; --i) {
    PointsEvent* cur_event = res_events[i];

    // If event is in between bridge (it was deleted), move one of bridge points.
    if (cur_event->pos.x > bridge_pt_1->pos.x &&
        cur_event->pos.x < bridge_pt_2->pos.x) {
      bridge_pt_1->next_point = bridge_pt_2->prev_point = cur_event;
      cur_event->prev_point = bridge_pt_1;
      cur_event->next_point = bridge_pt_2;
      if (cur_event->pos.x <= points[middle - 1].pos.x) {
        bridge_pt_1 = cur_event;
      } else {
        bridge_pt_2 = cur_event;
      }

      // Else, it was inserted and needs inverse.
    } else {
      cur_event->DeleteOrInsert();
      if (cur_event == bridge_pt_1) {
        bridge_pt_1 = bridge_pt_1->prev_point;
      }
      if (cur_event == bridge_pt_2) {
        bridge_pt_2 = bridge_pt_2->next_point;
      }
    }
  }
}


std::vector<PointsEvent*> GetEventsInHullBuilding(
    std::vector<LinkedPoint3>& points,
    const size_t begin,
    const size_t end) {
  if (end - begin == 1) {
    return std::vector<PointsEvent*>();
  }

  const size_t              middle = (begin + end) / 2;
  std::vector<PointsEvent*> res_events;
  std::vector<PointsEvent*> half_events[2] = {
    GetEventsInHullBuilding(points, begin, middle),
    GetEventsInHullBuilding(points, middle, end)
  };

  auto [bridge_pt_1, bridge_pt_2] = FindBridge(points, middle);

  size_t cur_event_1 = 0;
  size_t cur_event_2 = 0;
  for (double cur_time = -INF; ; ) {
    LinkedPoint3* cur_left_point  = nullptr;
    LinkedPoint3* cur_right_point = nullptr;
    double new_time[6]  = {};

    if (cur_event_1 < half_events[0].size()) {
      cur_left_point = half_events[0][cur_event_1];
      new_time[0] = GetTimeOfTurning(cur_left_point->prev_point,
                                     cur_left_point, cur_left_point->next_point);
    } else {
      new_time[0] = INF;
    }
    if (cur_event_2 < half_events[1].size()) {
      cur_right_point = half_events[1][cur_event_2];
      new_time[1] = GetTimeOfTurning(cur_right_point->prev_point,
                                     cur_right_point, cur_right_point->next_point);
    } else {
      new_time[1] = INF;
    }
    new_time[2] = GetTimeOfTurning(bridge_pt_1, bridge_pt_2, bridge_pt_2->next_point);
    new_time[3] = GetTimeOfTurning(bridge_pt_1, bridge_pt_2->prev_point, bridge_pt_2);
    new_time[4] = GetTimeOfTurning(bridge_pt_1->prev_point, bridge_pt_1, bridge_pt_2);
    new_time[5] = GetTimeOfTurning(bridge_pt_1, bridge_pt_1->next_point, bridge_pt_2);

    size_t min_time_ind = 7;
    double min_time     = INF;
    for (int i = 0; i < 6; i++) {
      if (new_time[i] > cur_time && new_time[i] < min_time) {
        min_time = new_time[i];
        min_time_ind = i;
      }
    }
    if (min_time_ind == 7) {
      break;
    }

    switch (min_time_ind) {
    case 0:
      if (cur_left_point->pos.x < bridge_pt_1->pos.x) {
        res_events.push_back(cur_left_point);
      }
      cur_left_point->DeleteOrInsert();
      cur_event_1++;
      break;
    case 1:
      if (cur_right_point->pos.x > bridge_pt_2->pos.x) {
        res_events.push_back(cur_right_point);
      }
      cur_right_point->DeleteOrInsert();
      cur_event_2++;
      break;
    case 2:
      res_events.push_back(bridge_pt_2);
      bridge_pt_2 = bridge_pt_2->next_point;
      break;
    case 3:
      bridge_pt_2 = bridge_pt_2->prev_point;
      res_events.push_back(bridge_pt_2);
      break;
    case 4:
      res_events.push_back(bridge_pt_1);
      bridge_pt_1 = bridge_pt_1->prev_point;
      break;
    case 5:
      bridge_pt_1 = bridge_pt_1->next_point;
      res_events.push_back(bridge_pt_1);
      break;
    default:
      break;
    }

    cur_time = min_time;
  }

  ReverseDeletionAndInsertion(points, middle, res_events, bridge_pt_1, bridge_pt_2);

  return res_events;
}


void BuildLowerOrUpperConvexHull(
    std::vector<LinkedPoint3>& points,
    Hull& res_convex_hull,
    bool lower_hull) {
  std::vector<PointsEvent*> events =
      GetEventsInHullBuilding(points, 0, points.size());
  for (auto event : events) {
    Face current = {event->prev_point->id, event->id, event->next_point->id};

    // If we need to delete (insert) a point, then points are ordered
    // clockwise relative to normal to exterior side, so make them counter-clockwise.
    if (lower_hull ? !event->DeleteOrInsert() : event->DeleteOrInsert()) {
      std::swap(current.point_id[0], current.point_id[1]);
    }
    res_convex_hull.faces.push_back(current);
  }
}


Hull BuildConvexHull(std::vector<LinkedPoint3> points) {
  const size_t points_count = points.size();
  Hull convex_hull = {};
  convex_hull.faces.reserve(points_count);

  // Sorting by x coordinate.
  std::sort(points.begin(), points.end());

  // Building lower convex hull, without flipping by z coordinate.
  BuildLowerOrUpperConvexHull(points, convex_hull, true);

  // Flip z coordinate to build upper hull.
  for (LinkedPoint3& p : points) {
    p.next_point = nullptr;
    p.prev_point = nullptr;
    p.pos.z = -p.pos.z;
  }
  BuildLowerOrUpperConvexHull(points, convex_hull, false);

  return convex_hull;
}


bool LinkedPoint3::DeleteOrInsert() {
  if (prev_point->next_point == this) {
    prev_point->next_point = next_point;
    next_point->prev_point = prev_point;
    return false;
  }
  prev_point->next_point = this;
  next_point->prev_point = this;
  return true;
}


Face::Face(std::initializer_list<size_t> list) {
  size_t i = 0;
  for (auto& elem : list) {
    point_id[i++] = elem;
  }
}


void Face::MakeMinIdPointFirst() {
  if (point_id[2] < point_id[0] && point_id[2] < point_id[1]) {
    std::swap(point_id[1], point_id[2]);
    std::swap(point_id[0], point_id[1]);
  } else if (point_id[1] < point_id[0] && point_id[1] < point_id[2]) {
    std::swap(point_id[0], point_id[1]);
    std::swap(point_id[1], point_id[2]);
  }
}


void Hull::SortFacesLexicographically() {
  for (auto& face : faces) {
    face.MakeMinIdPointFirst();
  }

  std::sort(faces.begin(), faces.end());
}


double GetCrossProductOfPoints(
    const LinkedPoint3* point_1,
    const LinkedPoint3* point_2,
    const LinkedPoint3* point_3) {
  if (point_1 == nullptr || point_2 == nullptr || point_3 == nullptr) {
    return INF;
  }
  return GetCrossProductXY(point_2->pos - point_1->pos, point_3->pos - point_2->pos);
}


double GetTimeOfTurning(
    const LinkedPoint3* point_1,
    const LinkedPoint3* point_2,
    const LinkedPoint3* point_3) {
  if (point_1 == nullptr || point_2 == nullptr || point_3 == nullptr) {
    return INF;
  }
  return GetCrossProductXZ(point_2->pos - point_1->pos, point_3->pos - point_2->pos) /
         GetCrossProductOfPoints(point_1, point_2, point_3);
}


bool operator<(const LinkedPoint3& p1, const LinkedPoint3& p2) {
  return p1.pos.x < p2.pos.x;
}


bool operator<(const Face& f1, const Face& f2) {
  if (f1.point_id[0] < f2.point_id[0]) {
    return true;
  } else if (f1.point_id[0] > f2.point_id[0]) {
    return false;
  }

  // f1[0] == f2[0];
  if (f1.point_id[1] < f2.point_id[1]) {
    return true;
  } else if (f1.point_id[1] > f2.point_id[1]) {
    return false;
  }

  // f1[1] == f2[1] && f1[0] == f2[0];
  return f1.point_id[2] < f2.point_id[2];
}


std::ostream& operator<<(std::ostream& stream, const Hull& hull) {
  stream << hull.faces.size() << std::endl;
  for (auto& face : hull.faces) {
    std::cout << 3                 << " " << 
                 face.point_id[0]  << " " <<
                 face.point_id[1]  << " " <<
                 face.point_id[2]  << std::endl;
  }

  return stream;
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


double GetCrossProductXY(const Vector3& vec_1, const Vector3& vec_2) {
  return vec_1.x * vec_2.y - vec_1.y * vec_2.x;
}


double GetCrossProductXZ(const Vector3& vec_1, const Vector3& vec_2) {
  return vec_1.x * vec_2.z - vec_1.z * vec_2.x;
}


void Vector3::RotateOX(double radians) {
  const double z_new = z * cos(radians) + y * sin(radians);
  const double y_new = -z * sin(radians) + y * cos(radians);
  z = z_new;
  y = y_new;
}


void Vector3::RotateOY(double radians) {
  const double x_new = x * cos(radians) + z * sin(radians);
  const double z_new = -x * sin(radians) + z * cos(radians);
  x = x_new;
  z = z_new;
}


void Vector3::RotateOZ(double radians) {
  const double x_new = x * cos(radians) + y * sin(radians);
  const double y_new = -x * sin(radians) + y * cos(radians);
  x = x_new;
  y = y_new;
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
