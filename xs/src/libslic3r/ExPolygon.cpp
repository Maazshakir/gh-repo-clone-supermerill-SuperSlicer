#include "BoundingBox.hpp"
#include "ExPolygon.hpp"
#include "Geometry.hpp"
#include "Polygon.hpp"
#include "Line.hpp"
#include "ClipperUtils.hpp"
#include "SVG.hpp"
#include "polypartition.h"
#include "poly2tri/poly2tri.h"
#include <algorithm>
#include <cassert>
#include <list>

namespace Slic3r {

ExPolygon::operator Points() const
{
    Points points;
    Polygons pp = *this;
    for (Polygons::const_iterator poly = pp.begin(); poly != pp.end(); ++poly) {
        for (Points::const_iterator point = poly->points.begin(); point != poly->points.end(); ++point)
            points.push_back(*point);
    }
    return points;
}

ExPolygon::operator Polygons() const
{
    return to_polygons(*this);
}

ExPolygon::operator Polylines() const
{
    return to_polylines(*this);
}

void
ExPolygon::scale(double factor)
{
    contour.scale(factor);
    for (Polygons::iterator it = holes.begin(); it != holes.end(); ++it) {
        (*it).scale(factor);
    }
}

void
ExPolygon::translate(double x, double y)
{
    contour.translate(x, y);
    for (Polygons::iterator it = holes.begin(); it != holes.end(); ++it) {
        (*it).translate(x, y);
    }
}

void
ExPolygon::rotate(double angle)
{
    contour.rotate(angle);
    for (Polygons::iterator it = holes.begin(); it != holes.end(); ++it) {
        (*it).rotate(angle);
    }
}

void
ExPolygon::rotate(double angle, const Point &center)
{
    contour.rotate(angle, center);
    for (Polygons::iterator it = holes.begin(); it != holes.end(); ++it) {
        (*it).rotate(angle, center);
    }
}

double
ExPolygon::area() const
{
    double a = this->contour.area();
    for (Polygons::const_iterator it = this->holes.begin(); it != this->holes.end(); ++it) {
        a -= -(*it).area();  // holes have negative area
    }
    return a;
}

bool
ExPolygon::is_valid() const
{
    if (!this->contour.is_valid() || !this->contour.is_counter_clockwise()) return false;
    for (Polygons::const_iterator it = this->holes.begin(); it != this->holes.end(); ++it) {
        if (!(*it).is_valid() || (*it).is_counter_clockwise()) return false;
    }
    return true;
}

bool
ExPolygon::contains(const Line &line) const
{
    return this->contains((Polyline)line);
}

bool
ExPolygon::contains(const Polyline &polyline) const
{
    return diff_pl((Polylines)polyline, *this).empty();
}

bool
ExPolygon::contains(const Polylines &polylines) const
{
    #if 0
    BoundingBox bbox = get_extents(polylines);
    bbox.merge(get_extents(*this));
    SVG svg(debug_out_path("ExPolygon_contains.svg"), bbox);
    svg.draw(*this);
    svg.draw_outline(*this);
    svg.draw(polylines, "blue");
    #endif
    Polylines pl_out = diff_pl(polylines, *this);
    #if 0
    svg.draw(pl_out, "red");
    #endif
    return pl_out.empty();
}

bool
ExPolygon::contains(const Point &point) const
{
    if (!this->contour.contains(point)) return false;
    for (Polygons::const_iterator it = this->holes.begin(); it != this->holes.end(); ++it) {
        if (it->contains(point)) return false;
    }
    return true;
}

// inclusive version of contains() that also checks whether point is on boundaries
bool
ExPolygon::contains_b(const Point &point) const
{
    return this->contains(point) || this->has_boundary_point(point);
}

bool
ExPolygon::has_boundary_point(const Point &point) const
{
    if (this->contour.has_boundary_point(point)) return true;
    for (Polygons::const_iterator h = this->holes.begin(); h != this->holes.end(); ++h) {
        if (h->has_boundary_point(point)) return true;
    }
    return false;
}

bool
ExPolygon::overlaps(const ExPolygon &other) const
{
    #if 0
    BoundingBox bbox = get_extents(other);
    bbox.merge(get_extents(*this));
    static int iRun = 0;
    SVG svg(debug_out_path("ExPolygon_overlaps-%d.svg", iRun ++), bbox);
    svg.draw(*this);
    svg.draw_outline(*this);
    svg.draw_outline(other, "blue");
    #endif
    Polylines pl_out = intersection_pl((Polylines)other, *this);
    #if 0
    svg.draw(pl_out, "red");
    #endif
    if (! pl_out.empty())
        return true; 
    return ! other.contour.points.empty() && this->contains_b(other.contour.points.front());
}

void ExPolygon::simplify_p(double tolerance, Polygons* polygons) const
{
    Polygons pp = this->simplify_p(tolerance);
    polygons->insert(polygons->end(), pp.begin(), pp.end());
}

Polygons ExPolygon::simplify_p(double tolerance) const
{
    Polygons pp;
    pp.reserve(this->holes.size() + 1);
    // contour
    {
        Polygon p = this->contour;
        p.points.push_back(p.points.front());
        p.points = MultiPoint::_douglas_peucker(p.points, tolerance);
        p.points.pop_back();
        pp.emplace_back(std::move(p));
    }
    // holes
    for (Polygon p : this->holes) {
        p.points.push_back(p.points.front());
        p.points = MultiPoint::_douglas_peucker(p.points, tolerance);
        p.points.pop_back();
        pp.emplace_back(std::move(p));
    }
    return simplify_polygons(pp);
}

ExPolygons ExPolygon::simplify(double tolerance) const
{
    return union_ex(this->simplify_p(tolerance));
}

void ExPolygon::simplify(double tolerance, ExPolygons* expolygons) const
{
    append(*expolygons, this->simplify(tolerance));
}

Point
ExPolygon::centroid() const{
    double cx = 0, cy = 0, double_area = 0;
    for (int i = 0; i < this->contour.points.size() - 1; i++) {
        const double mult = (this->contour.points[i].x * this->contour.points[i + 1].y) - (this->contour.points[i + 1].x * this->contour.points[i].y);
        cx += (this->contour.points[i].x + this->contour.points[i + 1].x) * mult;
        cy += (this->contour.points[i].y + this->contour.points[i + 1].y) * mult;
        double_area += mult;
    }
    double div = 1 / (3 * double_area);
    return Point((coord_t)(cx *div), (coord_t)(cy *div));
}

/// remove point that are at SCALED_EPSILON * 2 distance.
void
ExPolygon::remove_point_too_near(const coord_t smallest) {
    size_t id = 1;
    while (id < this->contour.points.size() - 1) {
        size_t newdist = min(this->contour.points[id].distance_to(this->contour.points[id - 1])
            , this->contour.points[id].distance_to(this->contour.points[id + 1]));
        if (newdist < smallest) {
            this->contour.points.erase(this->contour.points.begin() + id);
            newdist = this->contour.points[id].distance_to(this->contour.points[id - 1]);
        }
        //go to next one
        //if you removed a point, it check if the next one isn't too near from the previous one.
        // if not, it byepass it.
        if (newdist > smallest) {
            ++id;
        }
    }
}

//------- functions for medial_axis -----------

Point
interpolate(double percent, Polyline& poly) {
    const double pattern_length = poly.length();
    const double percent_epsilon = SCALED_EPSILON / pattern_length;
    double percent_lengthm1 = 0;
    double percent_length = 0;
    for (size_t idx_point = 1; idx_point < poly.points.size(); ++idx_point) {
        percent_lengthm1 = percent_length;
        percent_length += poly.points[idx_point - 1].distance_to(poly.points[idx_point]) / pattern_length;

        if (percent < percent_length + percent_epsilon) {
            //insert a new point before the position
            percent_length = (percent - percent_lengthm1) / (percent_length - percent_lengthm1);
            return poly.points[idx_point - 1].interpolate_to(percent_length, poly.points[idx_point]);
        }
    }
}



/// add points  from pattern to to_modify at the same % of the length
/// so not add if an other point is present at the correct position
void
add_point_same_percent(ThickPolyline* pattern, ThickPolyline* to_modify)
{
    const double to_modify_length = to_modify->length();
    const double percent_epsilon = SCALED_EPSILON / to_modify_length;
    const double pattern_length = pattern->length();

    double percent_length = 0;
    for (size_t idx_point = 1; idx_point < pattern->points.size() - 1; ++idx_point) {
        percent_length += pattern->points[idx_point-1].distance_to(pattern->points[idx_point]) / pattern_length;
        //find position 
        size_t idx_other = 1;
        double percent_length_other_before = 0;
        double percent_length_other = 0;
        while (idx_other < to_modify->points.size()) {
            percent_length_other_before = percent_length_other;
            percent_length_other += to_modify->points[idx_other-1].distance_to(to_modify->points[idx_other])
                / to_modify_length;
            if (percent_length_other > percent_length - percent_epsilon) {
                //if higher (we have gone over it)
                break;
            }
            ++idx_other;
        }
        if (percent_length_other > percent_length + percent_epsilon) {
            //insert a new point before the position
            double percent_dist = (percent_length - percent_length_other_before) / (percent_length_other - percent_length_other_before);
            coordf_t new_width = to_modify->width[idx_other - 1] * (1 - percent_dist)
                + to_modify->width[idx_other] * (percent_dist);
            to_modify->width.insert(to_modify->width.begin() + idx_other, new_width);
            to_modify->points.insert(to_modify->points.begin() + idx_other, 
                to_modify->points[idx_other - 1].interpolate_to(percent_dist, to_modify->points[idx_other]) );
        }
    }
}

/// find the nearest angle in the contour (or 2 nearest if it's difficult to choose) 
/// return 1 for an angle of 90° and 0 for an angle of 0° or 180°
double
get_coeff_from_angle_countour(Point &point, const ExPolygon &contour, coord_t min_dist_between_point)
{
    double nearestDist = point.distance_to(contour.contour.points.front());
    Point nearest = contour.contour.points.front();
    size_t id_nearest = 0;
    double nearDist = nearestDist;
    Point near = nearest;
    size_t id_near = 0;
    for (size_t id_point = 1; id_point < contour.contour.points.size(); ++id_point) {
        if (nearestDist > point.distance_to(contour.contour.points[id_point])) {
            //update near
            id_near = id_nearest;
            near = nearest;
            nearDist = nearestDist;
            //update nearest
            nearestDist = point.distance_to(contour.contour.points[id_point]);
            nearest = contour.contour.points[id_point];
            id_nearest = id_point;
        }
    }
    double angle = 0;
    size_t id_before = id_nearest == 0 ? contour.contour.points.size() - 1 : id_nearest - 1;
    Point point_before = id_nearest == 0 ? contour.contour.points.back() : contour.contour.points[id_nearest - 1];
    //Search one point far enough to be relevant
    while (nearest.distance_to(point_before) < min_dist_between_point) {
        point_before = id_before == 0 ? contour.contour.points.back() : contour.contour.points[id_before - 1];
        id_before = id_before == 0 ? contour.contour.points.size() - 1 : id_before - 1;
        //don't loop
        if (id_before == id_nearest) {
            id_before = id_nearest == 0 ? contour.contour.points.size() - 1 : id_nearest - 1;
            point_before = id_nearest == 0 ? contour.contour.points.back() : contour.contour.points[id_nearest - 1];
            break;
        }
    }
    size_t id_after = id_nearest == contour.contour.points.size() - 1 ? 0 : id_nearest + 1;
    Point point_after = id_nearest == contour.contour.points.size() - 1 ? contour.contour.points.front() : contour.contour.points[id_nearest + 1];
    //Search one point far enough to be relevant
    while (nearest.distance_to(point_after) < min_dist_between_point) {
        point_after = id_after == contour.contour.points.size() - 1 ? contour.contour.points.front() : contour.contour.points[id_after + 1];
        id_after = id_after == contour.contour.points.size() - 1 ? 0 : id_after + 1;
        //don't loop
        if (id_after == id_nearest) {
            id_after = id_nearest == contour.contour.points.size() - 1 ? 0 : id_nearest + 1;
            point_after = id_nearest == contour.contour.points.size() - 1 ? contour.contour.points.front() : contour.contour.points[id_nearest + 1];
            break;
        }
    }
    //compute angle
    angle = nearest.ccw_angle(point_before, point_after);
    if (angle >= PI) angle = 2*PI - angle;  // smaller angle
    //compute the diff from 90°
    angle = abs(angle - PI / 2);
    if (!near.coincides_with(nearest) && max(nearestDist, nearDist) + SCALED_EPSILON < nearest.distance_to(near)) {
        //not only nearest
        Point point_before = id_near == 0 ? contour.contour.points.back() : contour.contour.points[id_near - 1];
        Point point_after = id_near == contour.contour.points.size() - 1 ? contour.contour.points.front() : contour.contour.points[id_near + 1];
        double angle2 = min(nearest.ccw_angle(point_before, point_after), nearest.ccw_angle(point_after, point_before));
        angle2 = abs(angle - PI / 2);
        angle = (angle + angle2) / 2;
    }

    return 1-(angle/(PI/2));
}

void
fusion_curve(ExPolygon& polygon, ThickPolylines &pp, double max_width, double min_width) {

    //fusion Y with only 1 '0' value => the "0" branch "pull" the cross-point
    bool changes = false;
    for (size_t i = 0; i < pp.size(); ++i) {
        ThickPolyline& polyline = pp[i];
        std::cout
            << ", size: " << polyline.points.size()
            << ", length: " << unscale(polyline.length())
            << ", endf: " << polyline.endpoints.first
            << ", endb: " << polyline.endpoints.second
            << ", widthf: " << polyline.width.front()
            << ", widthb: " << polyline.width.back()
            << "\n";
        // only consider 2-point polyline with endpoint
        if (polyline.points.size() != 2) continue;
        if (polyline.endpoints.first) polyline.reverse();
        else if (!polyline.endpoints.second) continue;
        if (polyline.width.back() > 0) continue;

        //check my length is small
        coord_t length = (coord_t)polyline.length();
        if (length > max_width) continue;

        //only consider very shallow angle for contour
        if ( get_coeff_from_angle_countour(polyline.points.back(), polygon, SCALED_RESOLUTION) > 0.2) continue;

        // look if other end is a cross point with almost 90° angle
        double coeff_mainline = 0;
        // look if other end is a cross point with multiple other branch
        vector<size_t> crosspoint;
        for (size_t j = 0; j < pp.size(); ++j) {
            if (j == i) continue;
            ThickPolyline& other = pp[j];
            if (polyline.first_point().coincides_with(other.last_point())) {
                other.reverse();
                crosspoint.push_back(j);
                std::cout << " dot mainline test1: " << Line(polyline.points[0], polyline.points[1]).dot(Line(other.points[0], other.points[1]))
                    << ", p1=" << polyline.points[0].x << ":p2=" << polyline.points[1].x
                    << ", op1=" << other.points[0].x << ":op2=" << other.points[1].x
                    << ", p1=" << polyline.points[0].y << ":p2=" << polyline.points[1].y
                    << ", op1=" << other.points[0].y << ":op2=" << other.points[1].y
                    << "\n";
                coeff_mainline = max(coeff_mainline, 1 - abs(Line(polyline.points[0], polyline.points[1]).dot(Line(other.points[0], other.points[1]))));
            } else if (polyline.first_point().coincides_with(other.first_point())) {
                crosspoint.push_back(j);
                std::cout << " dot mainline test2: " << Line(polyline.points[0], polyline.points[1]).dot(Line(other.points[0], other.points[1])) << "\n";
                coeff_mainline = max(coeff_mainline, 1 - abs(Line(polyline.points[0], polyline.points[1]).dot(Line(other.points[0], other.points[1]))));
            }
        }
        //check if it's a line that we can pull
        //if (crosspoint.size() != 2) continue;
        std::cout << " coeff_mainline: " << coeff_mainline << "\n";
        if (coeff_mainline < 0.8) continue;

        // pull it a bit, depends on my size, the dot?, and the coeff at my 0-end (~12% for a square, almost 0 for a gentle curve)
        
        coord_t length_pull = polyline.length();
        length_pull *= 0.15 * get_coeff_from_angle_countour(polyline.points.back(), polygon, min(min_width, polyline.length() / 2));

        //compute dir
        Vectorf pull_direction(polyline.points[1].x - polyline.points[0].x, polyline.points[1].y - polyline.points[0].y);
        pull_direction = normalize(pull_direction);
        pull_direction.x *= length_pull;
        pull_direction.y *= length_pull;

        //pull the points
        Point &p1 = pp[crosspoint[0]].points[0];
        p1.x = p1.x + (coord_t)pull_direction.x;
        p1.y = p1.y + (coord_t)pull_direction.y;

        Point &p2 = pp[crosspoint[1]].points[0];
        p2.x = p2.x + (coord_t)pull_direction.x;
        p2.y = p2.y + (coord_t)pull_direction.y;
        

        std::cout << " remove curve: " << i << " / " << pp.size() << "\n";
        //delete the now unused polyline
        pp.erase(pp.begin() + i);
        --i;
        changes = true;
    }
    if (changes) {
        concatThickPolylines(pp);
        ///reorder, in case of change
        std::sort(pp.begin(), pp.end(), [](const ThickPolyline & a, const ThickPolyline & b) { return a.length() < b.length(); });
    }
}

void
fusion_corners(ExPolygon& polygon, ThickPolylines &pp, double max_width, double min_width) {

    //fusion Y with only 1 '0' value => the "0" branch "pull" the cross-point
    bool changes = false;
    for (size_t i = 0; i < pp.size(); ++i) {
        ThickPolyline& polyline = pp[i];
        // only consider polyline with 0-end
        if (polyline.points.size() != 2) continue;
        if (polyline.endpoints.first) polyline.reverse();
        else if (!polyline.endpoints.second) continue;
        if (polyline.width.back() > 0) continue;

        //check my length is small
        coord_t length = polyline.length();
        if (length > max_width) continue;

        // look if other end is a cross point with multiple other branch
        vector<size_t> crosspoint;
        for (size_t j = 0; j < pp.size(); ++j) {
            if (j == i) continue;
            ThickPolyline& other = pp[j];
            if (polyline.first_point().coincides_with(other.last_point())) {
                other.reverse();
                crosspoint.push_back(j);
            } else if (polyline.first_point().coincides_with(other.first_point())) {
                crosspoint.push_back(j);
            }
        }
        //check if it's a line that we can pull
        if (crosspoint.size() != 2) continue;

        // check if i am at the external side of a curve
        double angle1 = polyline.points[0].ccw_angle(polyline.points[1], pp[crosspoint[0]].points[1]);
        if (angle1 >= PI) angle1 = 2 * PI - angle1;  // smaller angle
        double angle2 = polyline.points[0].ccw_angle(polyline.points[1], pp[crosspoint[1]].points[1]);
        if (angle2 >= PI) angle2 = 2 * PI - angle2;  // smaller angle
        if (angle1 + angle2 < PI) continue;

        //check if is smaller or the other ones are not endpoits
        if (pp[crosspoint[0]].endpoints.second && length > pp[crosspoint[0]].length()) continue;
        if (pp[crosspoint[1]].endpoints.second && length > pp[crosspoint[1]].length()) continue;

        // if true, pull it a bit, depends on my size, the dot?, and the coeff at my 0-end (~12% for a square, almost 0 for a gentle curve)
        coord_t length_pull = polyline.length();
        length_pull *= 0.15 * get_coeff_from_angle_countour(polyline.points.back(), polygon, min(min_width, polyline.length() / 2));

        //compute dir
        Vectorf pull_direction(polyline.points[1].x - polyline.points[0].x, polyline.points[1].y - polyline.points[0].y);
        pull_direction = normalize(pull_direction);
        pull_direction.x *= length_pull;
        pull_direction.y *= length_pull;

        //pull the points
        Point &p1 = pp[crosspoint[0]].points[0];
        p1.x = p1.x + pull_direction.x;
        p1.y = p1.y + pull_direction.y;

        Point &p2 = pp[crosspoint[1]].points[0];
        p2.x = p2.x + pull_direction.x;
        p2.y = p2.y + pull_direction.y;

        //delete the now unused polyline
        pp.erase(pp.begin() + i);
        --i;
        changes = true;
    }
    if (changes) {
        concatThickPolylines(pp);
        ///reorder, in case of change
        std::sort(pp.begin(), pp.end(), [](const ThickPolyline & a, const ThickPolyline & b) { return a.length() < b.length(); });
    }
}


void
extends_line(ThickPolyline& polyline, const ExPolygon& inner_bound, const ExPolygon& outter_bounds, const ExPolygons& anchors, const coord_t max_width, const coord_t join_width) {

    // extend initial and final segments of each polyline if they're actual endpoints
    // We assign new endpoints to temporary variables because in case of a single-line
    // polyline, after we extend the start point it will be caught by the intersection()
    // call, so we keep the inner point until we perform the second intersection() as well
    std::cout << " ok, polyline.endpoints.second:" << polyline.endpoints.second << ", bounds?" << outter_bounds.has_boundary_point(polyline.points.back()) << "\n";
    if (polyline.endpoints.second && !outter_bounds.has_boundary_point(polyline.points.back())) {
        Line line(*(polyline.points.end() - 2), polyline.points.back());
        std::cout << " ok, line:" << unscale(line.a.x) << "->" << unscale(line.b.x) << " (x)\n";

        // prevent the line from touching on the other side, otherwise intersection() might return that solution
        if (polyline.points.size() == 2) line.a = line.midpoint();

        line.extend_end(max_width);
        Point new_back;
        if (inner_bound.contour.has_boundary_point(polyline.points.back())) {
            std::cout << "inner boundary ok\n";
            new_back = polyline.points.back();
        } else {
            (void)inner_bound.contour.first_intersection(line, &new_back);
            polyline.points.push_back(new_back);
            polyline.width.push_back(polyline.width.back());
            std::cout << "extends inner boundary to" << unscale(new_back.x) << "\n";
        }
        Point new_bound;
        (void)outter_bounds.contour.first_intersection(line, &new_bound);
        if (new_bound.coincides_with_epsilon(new_back)) {
            std::cout << "outter boundary ok\n";
            return;
        }
        //find anchor
        Point best_anchor;
        double shortest_dist = max_width;
        for (const ExPolygon& a : anchors) {
            std::cout << "try anchor\n";
            Point p_maybe_inside = a.centroid();
            std::cout << "total point " << unscale(p_maybe_inside.x) << "\n";
            double test_dist = new_bound.distance_to(p_maybe_inside) + new_back.distance_to(p_maybe_inside);
            std::cout << "test_dist= " << unscale(test_dist) << "? " << unscale(max_width)<<"\n";
            //if (test_dist < max_width / 2 && (test_dist < shortest_dist || shortest_dist < 0)) {
            if (test_dist < max_width && test_dist<shortest_dist) {
                std::cout << " best one! " << "\n";
                shortest_dist = test_dist;
                best_anchor = p_maybe_inside + new_bound;
                Line l2 = Line(line.a, p_maybe_inside);
                l2.extend_end(max_width);
                (void)outter_bounds.contour.first_intersection(l2, &new_bound);
            }
        }
        polyline.points.push_back(new_bound);
        polyline.width.push_back(join_width);
    }
}

void
fusion(const ExPolygon &bounds, double max_width, double min_width, ThickPolylines& pp, ExPolygon& polygon)
{
    // Aligned fusion: Fusion the bits at the end of lines by "increasing thickness"
    // For that, we have to find other lines,
    // and with a next point no more distant than the max width.
    // Then, we can merge the bit from the first point to the second by following the mean.
    //
    bool changes = true;
    while (changes) {
        changes = false;
        for (size_t i = 0; i < pp.size(); ++i) {
            ThickPolyline& polyline = pp[i];

            //simple check to see if i can be fusionned
            if (!polyline.endpoints.first && !polyline.endpoints.second) continue;


            ThickPolyline* best_candidate = nullptr;
            float best_dot = -1;
            size_t best_idx = 0;
            double dot_poly_branch = 0;
            double dot_candidate_branch = 0;

            // find another polyline starting here
            for (size_t j = i + 1; j < pp.size(); ++j) {
                ThickPolyline& other = pp[j];
                if (polyline.last_point().coincides_with(other.last_point())) {
                    polyline.reverse();
                    other.reverse();
                } else if (polyline.first_point().coincides_with(other.last_point())) {
                    other.reverse();
                } else if (polyline.first_point().coincides_with(other.first_point())) {
                } else if (polyline.last_point().coincides_with(other.first_point())) {
                    polyline.reverse();
                } else {
                    continue;
                }
                //std::cout << " try : " << i << ":" << j << " : " << 
                //    (polyline.points.size() < 2 && other.points.size() < 2) <<
                //    (!polyline.endpoints.second || !other.endpoints.second) <<
                //    ((polyline.points.back().distance_to(other.points.back())
                //    + (polyline.width.back() + other.width.back()) / 4)
                //    > max_width*1.05) <<
                //    (abs(polyline.length() - other.length()) > max_width / 2) << "\n";

                //// mergeable tests
                if (polyline.points.size() < 2 && other.points.size() < 2) continue;
                if (!polyline.endpoints.second || !other.endpoints.second) continue;
                // test if the new width will not be too big if a fusion occur
                //note that this isn't the real calcul. It's just to avoid merging lines too far apart.
                if (
                    ((polyline.points.back().distance_to(other.points.back())
                    + (polyline.width.back() + other.width.back()) / 4)
                > max_width*1.05))
                continue;
                // test if the lines are not too different in length.
                if (abs(polyline.length() - other.length()) > max_width) continue;


                //convex test
                Point garbage_point;
                if (other.length() > min_width &&
                    bounds.contour.intersection(Line(interpolate(0.9, polyline), interpolate(0.9, other)), &garbage_point)
                    ) continue;

                //compute angle to see if it's better than previous ones (straighter = better).
                const float other_dot = polyline.lines().front().dot(other.lines().front());

                // Get the branch/line in wich we may merge, if possible
                // with that, we can decide what is important, and how we can merge that.
                // angle_poly - angle_candi =90° => one is useless
                // both angle are equal => both are useful with same strength
                // ex: Y => | both are useful to crete a nice line
                // ex2: TTTTT => -----  these 90° useless lines should be discarded
                bool find_main_branch = false;
                size_t biggest_main_branch_id = 0;
                coord_t biggest_main_branch_length = 0;
                for (size_t k = 0; k < pp.size(); ++k) {
                    //std::cout << "try to find main : " << k << " ? " << i << " " << j << " ";
                    if (k == i | k == j) continue;
                    ThickPolyline& main = pp[k];
                    if (polyline.first_point().coincides_with(main.last_point())) {
                        main.reverse();
                        if (!main.endpoints.second)
                            find_main_branch = true;
                        else if (biggest_main_branch_length < main.length()) {
                            biggest_main_branch_id = k;
                            biggest_main_branch_length = main.length();
                        }
                    } else if (polyline.first_point().coincides_with(main.first_point())) {
                        if (!main.endpoints.second)
                            find_main_branch = true;
                        else if (biggest_main_branch_length < main.length()) {
                            biggest_main_branch_id = k;
                            biggest_main_branch_length = main.length();
                        }
                    }
                    if (find_main_branch) {
                        //use this variable to store the good index and break to compute it
                        biggest_main_branch_id = k;
                        break;
                    }
                }
                if (!find_main_branch && biggest_main_branch_length == 0) {
                    // nothing -> it's impossible!
                    dot_poly_branch = 0.707;
                    dot_candidate_branch = 0.707;
                    //std::cout << "no main branch... impossible!!\n";
                } else if (!find_main_branch &&
                    (pp[biggest_main_branch_id].length() < polyline.length() || pp[biggest_main_branch_id].length() < other.length())) {
                    //the main branch should have no endpoint or be bigger!
                    //here, it have an endpoint, and is not the biggest -> bad!
                    continue;
                } else {
                    //compute the dot (biggest_main_branch_id)
                    Pointf v_poly(polyline.lines().front().vector().x, polyline.lines().front().vector().y);
                    v_poly.scale(1 / std::sqrt(v_poly.x*v_poly.x + v_poly.y*v_poly.y));
                    Pointf v_candid(other.lines().front().vector().x, other.lines().front().vector().y);
                    v_candid.scale(1 / std::sqrt(v_candid.x*v_candid.x + v_candid.y*v_candid.y));
                    Pointf v_branch(-pp[biggest_main_branch_id].lines().front().vector().x, -pp[biggest_main_branch_id].lines().front().vector().y);
                    v_branch.scale(1 / std::sqrt(v_branch.x*v_branch.x + v_branch.y*v_branch.y));
                    dot_poly_branch = v_poly.x*v_branch.x + v_poly.y*v_branch.y;
                    dot_candidate_branch = v_candid.x*v_branch.x + v_candid.y*v_branch.y;

                    if (dot_poly_branch < 0) dot_poly_branch = 0;
                    if (dot_candidate_branch < 0) dot_candidate_branch = 0;
                }
                //test if it's useful to merge or not
                //ie, don't merge  'T' but ok for 'Y', merge only lines of not disproportionate different length (ratio max: 4)
                if (dot_poly_branch < 0.1 || dot_candidate_branch < 0.1 ||
                    (polyline.length()>other.length() ? polyline.length() / other.length() : other.length() / polyline.length()) > 4) {
                    continue;
                }
                if (other_dot > best_dot) {
                    best_candidate = &other;
                    best_idx = j;
                    best_dot = other_dot;
                }
            }
            if (best_candidate != nullptr) {
                // delete very near points
                polyline.remove_point_too_near();
                best_candidate->remove_point_too_near();

                // add point at the same pos than the other line to have a nicer fusion
                add_point_same_percent(&polyline, best_candidate);
                add_point_same_percent(best_candidate, &polyline);

                //get the angle of the nearest points of the contour to see : _| (good) \_ (average) __(bad)
                //sqrt because the result are nicer this way: don't over-penalize /_ angles
                //TODO: try if we can achieve a better result if we use a different algo if the angle is <90°
                const double coeff_angle_poly = (get_coeff_from_angle_countour(polyline.points.back(), polygon, min(min_width, polyline.length() / 2)));
                const double coeff_angle_candi = (get_coeff_from_angle_countour(best_candidate->points.back(), polygon, min(min_width, best_candidate->length() / 2)));

                //this will encourage to follow the curve, a little, because it's shorter near the center
                //without that, it tends to go to the outter rim.
                double weight_poly = 1 - polyline.length() / max(polyline.length(), best_candidate->length());
                double weight_candi = 1 - best_candidate->length() / max(polyline.length(), best_candidate->length());
                weight_poly *= 0;
                weight_candi *= 0;
                weight_poly += 1;
                weight_candi += 1;
                weight_poly = 0.5 + coeff_angle_poly;
                weight_candi = 0.5 + coeff_angle_candi;
                const double coeff_poly = (dot_poly_branch * weight_poly) / (dot_poly_branch * weight_poly + dot_candidate_branch * weight_candi);
                const double coeff_candi = 1.0 - coeff_poly;
                //iterate the points
                // as voronoi should create symetric thing, we can iterate synchonously
                size_t idx_point = 1;
                while (idx_point < min(polyline.points.size(), best_candidate->points.size())) {
                    //fusion
                    polyline.points[idx_point].x = polyline.points[idx_point].x * coeff_poly + best_candidate->points[idx_point].x * coeff_candi;
                    polyline.points[idx_point].y = polyline.points[idx_point].y * coeff_poly + best_candidate->points[idx_point].y * coeff_candi;

                    // The width decrease with distance from the centerline.
                    // This formula is what works the best, even if it's not perfect (created empirically).  0->3% error on a gap fill on some tests.
                    //If someone find  an other formula based on the properties of the voronoi algorithm used here, and it works better, please use it.
                    //or maybe just use the distance to nearest edge in bounds...
                    double value_from_current_width = 0.5*polyline.width[idx_point] * dot_poly_branch / max(dot_poly_branch, dot_candidate_branch);
                    value_from_current_width += 0.5*best_candidate->width[idx_point] * dot_candidate_branch / max(dot_poly_branch, dot_candidate_branch);
                    double value_from_dist = 2 * polyline.points[idx_point].distance_to(best_candidate->points[idx_point]);
                    value_from_dist *= sqrt(min(dot_poly_branch, dot_candidate_branch) / max(dot_poly_branch, dot_candidate_branch));
                    polyline.width[idx_point] = value_from_current_width + value_from_dist;
                    //failsafe
                    if (polyline.width[idx_point] > max_width) polyline.width[idx_point] = max_width;

                    ++idx_point;
                }
                if (idx_point < best_candidate->points.size()) {
                    if (idx_point + 1 < best_candidate->points.size()) {
                        //create a new polyline
                        pp.emplace_back();
                        pp.back().endpoints.first = true;
                        pp.back().endpoints.second = best_candidate->endpoints.second;
                        for (size_t idx_point_new_line = idx_point; idx_point_new_line < best_candidate->points.size(); ++idx_point_new_line) {
                            pp.back().points.push_back(best_candidate->points[idx_point_new_line]);
                            pp.back().width.push_back(best_candidate->width[idx_point_new_line]);
                        }
                    } else {
                        //Add last point
                        polyline.points.push_back(best_candidate->points[idx_point]);
                        polyline.width.push_back(best_candidate->width[idx_point]);
                        //select if an end opccur
                        polyline.endpoints.second &= best_candidate->endpoints.second;
                    }

                } else {
                    //select if an end opccur
                    polyline.endpoints.second &= best_candidate->endpoints.second;
                }

                //remove points that are the same or too close each other, ie simplify
                for (size_t idx_point = 1; idx_point < polyline.points.size(); ++idx_point) {
                    if (polyline.points[idx_point - 1].distance_to(polyline.points[idx_point]) < SCALED_EPSILON) {
                        if (idx_point < polyline.points.size() - 1) {
                            polyline.points.erase(polyline.points.begin() + idx_point);
                            polyline.width.erase(polyline.width.begin() + idx_point);
                        } else {
                            polyline.points.erase(polyline.points.begin() + idx_point - 1);
                            polyline.width.erase(polyline.width.begin() + idx_point - 1);
                        }
                        --idx_point;
                    }
                }
                //remove points that are outside of the geometry
                for (size_t idx_point = 0; idx_point < polyline.points.size(); ++idx_point) {
                    if (!bounds.contains_b(polyline.points[idx_point])) {
                        polyline.points.erase(polyline.points.begin() + idx_point);
                        polyline.width.erase(polyline.width.begin() + idx_point);
                        --idx_point;
                    }
                }
                if (polyline.points.size() < 2) {
                    //remove self
                    pp.erase(pp.begin() + i);
                    --i;
                    --best_idx;
                }

                pp.erase(pp.begin() + best_idx);
                changes = true;
                break;
            }
        }
        if (changes) {
            concatThickPolylines(pp);
            ///reorder, in case of change
            std::sort(pp.begin(), pp.end(), [](const ThickPolyline & a, const ThickPolyline & b) { return a.length() < b.length(); });
        }
    }
}

void
cut_too_thin_polylines(double max_width, double min_width, ThickPolylines& pp)
{
    bool changes = false;
    for (size_t i = 0; i < pp.size(); ++i) {
        ThickPolyline& polyline = pp[i];
        // remove bits with too small extrusion
        while (polyline.points.size() > 1 && polyline.width.front() < min_width && polyline.endpoints.first) {
            //try to split if possible
            if (polyline.width[1] > min_width) {
                double percent_can_keep = (min_width - polyline.width[0]) / (polyline.width[1] - polyline.width[0]);
                if (polyline.points.front().distance_to(polyline.points[1])* (1 - percent_can_keep) > max_width / 2) {
                    //Can split => move the first point and assign a new weight.
                    //the update of endpoints wil be performed in concatThickPolylines
                    polyline.points.front().x = polyline.points.front().x +
                        (coord_t)((polyline.points[1].x - polyline.points.front().x) * percent_can_keep);
                    polyline.points.front().y = polyline.points.front().y +
                        (coord_t)((polyline.points[1].y - polyline.points.front().y) * percent_can_keep);
                    polyline.width.front() = min_width;
                    changes = true;
                    break;
                }
            }
            polyline.points.erase(polyline.points.begin());
            polyline.width.erase(polyline.width.begin());
            changes = true;
        }
        while (polyline.points.size() > 1 && polyline.width.back() < min_width && polyline.endpoints.second) {
            //try to split if possible
            if (polyline.width[polyline.points.size() - 2] > min_width) {
                double percent_can_keep = (min_width - polyline.width.back()) / (polyline.width[polyline.points.size() - 2] - polyline.width.back());
                if (polyline.points.back().distance_to(polyline.points[polyline.points.size() - 2]) * (1 - percent_can_keep) > max_width / 2) {
                    //Can split => move the first point and assign a new weight.
                    //the update of endpoints wil be performed in concatThickPolylines
                    polyline.points.back().x = polyline.points.back().x +
                        (coord_t)((polyline.points[polyline.points.size() - 2].x - polyline.points.back().x) * percent_can_keep);
                    polyline.points.back().y = polyline.points.back().y +
                        (coord_t)((polyline.points[polyline.points.size() - 2].y - polyline.points.back().y) * percent_can_keep);
                    polyline.width.back() = min_width;
                    changes = true;
                    break;
                }
            }
            polyline.points.erase(polyline.points.end() - 1);
            polyline.width.erase(polyline.width.end() - 1);
            changes = true;
        }
        if (polyline.points.size() < 2) {
            //remove self if too small
            pp.erase(pp.begin() + i);
            --i;
        }
    }
    if (changes) concatThickPolylines(pp);

}

int id = 0;
void
ExPolygon::medial_axis(const ExPolygon &bounds, double max_width, double min_width, ThickPolylines* polylines, double height) const
{
    id++;
    std::cout << " ===== " << id << " ========\n";

    {
        stringstream stri;
        stri << "medial_axis_0.0_base_" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(*this);
        svg.Close();
    }

    //simplify the boundary between us and the bounds.
    //int firstIdx = 0;
    //while (firstIdx < contour.points.size() && bounds.contour.contains(contour.points[firstIdx])) firstIdx++;
    ExPolygon simplified_poly = *this;
    if (this != &bounds) {
        bool need_intersect = false;
        for (size_t i = 0; i < simplified_poly.contour.points.size(); i++) {
            Point &p_check = simplified_poly.contour.points[i];
            std::cout << "check point " << p_check.x << " : " << p_check.y << "\n";
            //bool find = false;
            //for (size_t j = 0; j<bounds.contour.points.size(); j++) {
            //    const Point &p_iter = bounds.contour.points[j];
            //    if (abs(p_iter.x - p_check.x)< SCALED_EPSILON
            //        && abs(p_iter.y - p_check.y)< SCALED_EPSILON) {
            //        std::cout << "find point " << "\n";
            //        find = true;
            //        break;
            //    }
            //    size_t next = (j + 1) % simplified_poly.contour.points.size();
            //    double angle = p_check.ccw_angle(p_iter, bounds.contour.points[next]);
            //    if (angle > PI) angle = PI * 2 - angle;
            //    if (angle > PI*0.99) {
            //        std::cout << "find point (almost inside a line) " << angle / PI << "\n";
            //        find = true;
            //        break;
            //    }
            //}
            //if (!find) {
            if (!bounds.has_boundary_point(p_check)) {
                //check if we put it at a bound point instead of delete it
                size_t prev_i = i == 0 ? simplified_poly.contour.points.size() - 1 : (i - 1);
                size_t next_i = i == simplified_poly.contour.points.size() - 1 ? 0 : (i + 1);
                const Point* closest = bounds.contour.closest_point(p_check);
                if (closest != nullptr && closest->distance_to(p_check) + SCALED_EPSILON
                    < min(p_check.distance_to(simplified_poly.contour.points[prev_i]), p_check.distance_to(simplified_poly.contour.points[next_i])) / 2) {
                    p_check.x = closest->x;
                    p_check.y = closest->y;
                    std::cout << " move to " << closest->x << " : " << closest->x << "\n";
                    need_intersect = true;
                } else {
                    std::cout << " erase " << "\n";
                    simplified_poly.contour.points.erase(simplified_poly.contour.points.begin() + i);
                    i--;
                }
            }
        }
        if (need_intersect) {
            ExPolygons simplified_polys = intersection_ex(simplified_poly, bounds);
            if (simplified_polys.size() == 1) {
                simplified_poly = simplified_polys[0];
            } else {
                simplified_poly = *this;
            }
        }
    }
    std::cout << "SCALED_RESOLUTION=" << unscale(SCALED_RESOLUTION) << ", /10=" << unscale(SCALED_RESOLUTION / 10)
        << ", SCALED_EPSILON=" << unscale(SCALED_EPSILON) << ", x5=" << unscale(SCALED_EPSILON * 10)
        << "\n";
    simplified_poly.remove_point_too_near(min(SCALED_RESOLUTION / 20, SCALED_EPSILON * 3));

    {
        stringstream stri;
        stri << "medial_axis_0.1_simplified_" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.Close();
    }


    // init helper object
    //Slic3r::Geometry::MedialAxis ma(max_width, min_width, &bounds);
    //ma.lines = bounds.lines();
    Slic3r::Geometry::MedialAxis ma(max_width, min_width, &simplified_poly);
    ma.lines = simplified_poly.lines();
    for (size_t i = 0; i < ma.lines.size(); i++) {
        Point bounds_p1 = ma.lines[i].a.projection_onto(bounds.contour);
        Point bounds_p2 = ma.lines[i].b.projection_onto(bounds.contour);
        //check if b1 and b2 are adjacent
        long idx_b1 = -20;
        long idx_b2 = -2;
        for (size_t j = 0; j < bounds.contour.points.size(); j++) {
            if (bounds.contour.points[j].coincides_with_epsilon(bounds_p1)) {
                idx_b1 = j;
            }
            if (bounds.contour.points[j].coincides_with_epsilon(bounds_p2)) {
                idx_b2 = j;
            }
            if (idx_b1 >= 0 && idx_b2 >= 0) break;
        }
        if (abs((idx_b1 - idx_b2)) > 1 && !(idx_b1 == 0 && idx_b2 == bounds.contour.points.size() - 1) && !(idx_b2 == 0 && idx_b1 == bounds.contour.points.size() - 1)) {
            std::cout << "del line " << idx_b1 << ", " << idx_b2 << " / " << bounds.contour.points.size()<< "  | " << i << " / " << ma.lines.size() << "\n";
            //not adjacent! remove it!
            ma.lines.erase(ma.lines.begin() + i);
            i--;
        }
    }
    
    // compute the Voronoi diagram and extract medial axis polylines
    ThickPolylines pp;
    ma.build(&pp);

    {
        stringstream stri;
        stri << "medial_axis_0.2_voronoi_" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(ma.lines);
        svg.draw(pp);
        svg.Close();
    }

    //{
    //    stringstream stri;
    //    stri << "medial_axis" << id << ".svg";
    //    SVG svg(stri.str());
    //    svg.draw(bounds);
    //    svg.draw(*this);
    //    svg.draw(pp);
    //    svg.Close();
    //}
    
    /* Find the maximum width returned; we're going to use this for validating and 
       filtering the output segments. */
    double max_w = 0;
    for (ThickPolylines::const_iterator it = pp.begin(); it != pp.end(); ++it)
        max_w = fmaxf(max_w, *std::max_element(it->width.begin(), it->width.end()));

    //reoder pp by length (ascending) It's really important to do that to avoid building the line from the width insteand of the length
    std::sort(pp.begin(), pp.end(), [](const ThickPolyline & a, const ThickPolyline & b) { return a.length() < b.length(); });
    {
        stringstream stri;
        stri << "medial_axis_1_concat" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.draw(pp);
        svg.Close();
    }
    for (ThickPolyline& p : pp) {
        std::cout << "polyline thick : ";
        for (double width : p.width) {
            std::cout << ", " << unscale(width);
        }
        std::cout << "\n";
    }


    fusion_curve(simplified_poly, pp, max_width, min_width);
    {
        stringstream stri;
        stri << "medial_axis_2_uncurved" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.draw(pp);
        svg.Close();
    }
    concatThickPolylines(pp);
    fusion(bounds, max_width, min_width, pp, simplified_poly);
    std::cout << "0\n";

    for (ThickPolyline& p : pp) {
        std::cout << "FUSIONED polyline thick : ";
        for (double width : p.width) {
            std::cout << ", " << unscale(width);
        }
        std::cout << "\n";
    }
    {
        stringstream stri;
        stri << "medial_axis_3fusioned" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.draw(pp);
        svg.Close();
    }
    fusion_corners(simplified_poly, pp, max_width, min_width);
    std::cout << "0\n";

    {
        stringstream stri;
        stri << "medial_axis_4fusioned_croner" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.draw(pp);
        svg.Close();
    }

    // remove too small extrusion at start & end of polylines
   

    // Loop through all returned polylines in order to extend their endpoints to the 
    //   expolygon boundaries
    ExPolygons anchors = diff_ex(bounds, *this, true);
    for (size_t i = 0; i < pp.size(); ++i) {
        ThickPolyline& polyline = pp[i];
        std::cout << "extends line1  " << i<<"\n";
        extends_line(polyline, simplified_poly, bounds, anchors, max_width, min_width);
        polyline.reverse();
        std::cout << "extends line2  " << i << "\n";
        extends_line(polyline, simplified_poly, bounds, anchors, max_width, min_width);
        std::cout << "extends line end  " << i << "\n";
    }

    {
        stringstream stri;
        stri << "medial_axis_5_extended" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.draw(pp);
        svg.Close();
    }

    cut_too_thin_polylines(max_width, min_width, pp);

    // concatenate, but even where multiple thickpolyline join, to create nice long strait polylines
    /*  If we removed any short polylines we now try to connect consecutive polylines
        in order to allow loop detection. Note that this algorithm is greedier than 
        MedialAxis::process_edge_neighbors() as it will connect random pairs of 
        polylines even when more than two start from the same point. This has no 
        drawbacks since we optimize later using nearest-neighbor which would do the 
        same, but should we use a more sophisticated optimization algorithm we should
        not connect polylines when more than two meet. 
        Optimisation of the old algorithm : now we select the most "strait line" choice 
        when we merge with an other line at a point with more than two meet.
        */
    bool changes = false;
    for (size_t i = 0; i < pp.size(); ++i) {
        ThickPolyline& polyline = pp[i];
        if (polyline.endpoints.first && polyline.endpoints.second) continue; // optimization
            
        ThickPolyline* best_candidate = nullptr;
        float best_dot = -1;
        size_t best_idx = 0;

        // find another polyline starting here
        for (size_t j = i + 1; j < pp.size(); ++j) {
            ThickPolyline& other = pp[j];
            if (polyline.last_point().coincides_with(other.last_point())) {
                other.reverse();
            } else if (polyline.first_point().coincides_with(other.last_point())) {
                polyline.reverse();
                other.reverse();
            } else if (polyline.first_point().coincides_with(other.first_point())) {
                polyline.reverse();
            } else if (!polyline.last_point().coincides_with(other.first_point())) {
                continue;
            }

            Pointf v_poly(polyline.lines().back().vector().x, polyline.lines().back().vector().y);
            v_poly.scale(1 / std::sqrt(v_poly.x*v_poly.x + v_poly.y*v_poly.y));
            Pointf v_other(other.lines().front().vector().x, other.lines().front().vector().y);
            v_other.scale(1 / std::sqrt(v_other.x*v_other.x + v_other.y*v_other.y));
            float other_dot = v_poly.x*v_other.x + v_poly.y*v_other.y;
            if (other_dot > best_dot) {
                best_candidate = &other;
                best_idx = j;
                best_dot = other_dot;
            }
        }
        if (best_candidate != nullptr) {

            polyline.points.insert(polyline.points.end(), best_candidate->points.begin() + 1, best_candidate->points.end());
            polyline.width.insert(polyline.width.end(), best_candidate->width.begin() + 1, best_candidate->width.end());
            polyline.endpoints.second = best_candidate->endpoints.second;
            assert(polyline.width.size() == polyline.points.size());
            changes = true;
            pp.erase(pp.begin() + best_idx);
        }
    }
    if (changes) concatThickPolylines(pp);

    //remove too thin polylines points (inside a polyline : split it)
    for (size_t i = 0; i < pp.size(); ++i) {
        ThickPolyline& polyline = pp[i];

        // remove bits with too small extrusion
        size_t idx_point = 0;
        while (idx_point<polyline.points.size()) {
            if (polyline.width[idx_point] < min_width) {
                if (idx_point == 0) {
                    //too thin at start
                    polyline.points.erase(polyline.points.begin());
                    polyline.width.erase(polyline.width.begin());
                    idx_point = 0;
                } else if (idx_point == 1) {
                    //too thin at start
                    polyline.points.erase(polyline.points.begin());
                    polyline.width.erase(polyline.width.begin());
                    polyline.points.erase(polyline.points.begin());
                    polyline.width.erase(polyline.width.begin());
                    idx_point = 0;
                } else if (idx_point == polyline.points.size() - 2) {
                    //too thin at (near) end
                    polyline.points.erase(polyline.points.end() - 1);
                    polyline.width.erase(polyline.width.end() - 1);
                    polyline.points.erase(polyline.points.end() - 1);
                    polyline.width.erase(polyline.width.end() - 1);
                } else if (idx_point == polyline.points.size() - 1) {
                    //too thin at end
                    polyline.points.erase(polyline.points.end() - 1);
                    polyline.width.erase(polyline.width.end() - 1);
                } else {
                    //too thin in middle : split
                    pp.emplace_back();
                    ThickPolyline &newone = pp.back();
                    newone.points.insert(newone.points.begin(), polyline.points.begin() + idx_point + 1, polyline.points.end());
                    newone.width.insert(newone.width.begin(), polyline.width.begin() + idx_point + 1, polyline.width.end());
                    polyline.points.erase(polyline.points.begin() + idx_point, polyline.points.end());
                    polyline.width.erase(polyline.width.begin() + idx_point, polyline.width.end());
                }
            } else idx_point++;

            if (polyline.points.size() < 2) {
                //remove self if too small
                pp.erase(pp.begin() + i);
                --i;
                break;
            }
        }
    }

    //remove too short polyline
    changes = true;
    while (changes) {
        changes = false;
        
        double shortest_size = max_w * 2;
        size_t shortest_idx = -1;
        for (size_t i = 0; i < pp.size(); ++i) {
            ThickPolyline& polyline = pp[i];
            // Remove the shortest polylines : polyline that are shorter than wider
            // (we can't do this check before endpoints extension and clipping because we don't
            // know how long will the endpoints be extended since it depends on polygon thickness
            // which is variable - extension will be <= max_width/2 on each side) 
            if ((polyline.endpoints.first || polyline.endpoints.second)
                && polyline.length() < max_width / 2) {
               if (shortest_size > polyline.length()) {
                    shortest_size = polyline.length();
                    shortest_idx = i;
                }

            }
        }
        if (shortest_idx >= 0 && shortest_idx < pp.size()) {
            pp.erase(pp.begin() + shortest_idx);
            changes = true;
        }
        if (changes) concatThickPolylines(pp);
    }

    //TODO: reduce the flow at the intersection ( + ) points ?

    //ensure the volume extruded is correct for what we have been asked
    // => don't over-extrude
    double surface = 0;
    double volume = 0;
    for (ThickPolyline& polyline : pp) {
        for (ThickLine l : polyline.thicklines()) {
            surface += l.length() * (l.a_width + l.b_width) / 2;
            double width_mean = (l.a_width + l.b_width) / 2;
            volume += height * (width_mean - height * (1. - 0.25 * PI)) * l.length();
        }
    }

    // compute bounds volume
    double boundsVolume = 0;
    boundsVolume += height*bounds.area();
    // add external "perimeter gap"
    double perimeterRoundGap = bounds.contour.length() * height * (1 - 0.25*PI) * 0.5;
    // add holes "perimeter gaps"
    double holesGaps = 0;
    for (auto hole = bounds.holes.begin(); hole != bounds.holes.end(); ++hole) {
        holesGaps += hole->length() * height * (1 - 0.25*PI) * 0.5;
    }
    boundsVolume += perimeterRoundGap + holesGaps;
    
    if (boundsVolume < volume) {
        //reduce width
        double reduce_by = boundsVolume / volume;
        for (ThickPolyline& polyline : pp) {
            for (ThickLine l : polyline.thicklines()) {
                l.a_width *= reduce_by;
                l.b_width *= reduce_by;
            }
        }
    }
    polylines->insert(polylines->end(), pp.begin(), pp.end());
    {
        stringstream stri;
        stri << "medial_axis_9_end" << id << ".svg";
        SVG svg(stri.str());
        svg.draw(bounds);
        svg.draw(simplified_poly);
        svg.draw(pp);
        svg.Close();
    }
    for (ThickPolyline& p : pp) {
        std::cout << "END polyline thick : ";
        for (double width : p.width) {
            std::cout << ", " << unscale(width);
        }
        std::cout << "\n";
    }
    std::cout << " END ===== " << id << " ========\n";
}

void
ExPolygon::medial_axis(double max_width, double min_width, Polylines* polylines) const
{
    ThickPolylines tp;
    this->medial_axis(*this, max_width, min_width, &tp, max_width/2.0);
    polylines->insert(polylines->end(), tp.begin(), tp.end());
}

void
ExPolygon::get_trapezoids(Polygons* polygons) const
{
    ExPolygons expp;
    expp.push_back(*this);
    boost::polygon::get_trapezoids(*polygons, expp);
}

void
ExPolygon::get_trapezoids(Polygons* polygons, double angle) const
{
    ExPolygon clone = *this;
    clone.rotate(PI/2 - angle, Point(0,0));
    clone.get_trapezoids(polygons);
    for (Polygons::iterator polygon = polygons->begin(); polygon != polygons->end(); ++polygon)
        polygon->rotate(-(PI/2 - angle), Point(0,0));
}

// This algorithm may return more trapezoids than necessary
// (i.e. it may break a single trapezoid in several because
// other parts of the object have x coordinates in the middle)
void
ExPolygon::get_trapezoids2(Polygons* polygons) const
{
    // get all points of this ExPolygon
    Points pp = *this;
    
    // build our bounding box
    BoundingBox bb(pp);
    
    // get all x coordinates
    std::vector<coord_t> xx;
    xx.reserve(pp.size());
    for (Points::const_iterator p = pp.begin(); p != pp.end(); ++p)
        xx.push_back(p->x);
    std::sort(xx.begin(), xx.end());
    
    // find trapezoids by looping from first to next-to-last coordinate
    for (std::vector<coord_t>::const_iterator x = xx.begin(); x != xx.end()-1; ++x) {
        coord_t next_x = *(x + 1);
        if (*x == next_x) continue;
        
        // build rectangle
        Polygon poly;
        poly.points.resize(4);
        poly[0].x = *x;
        poly[0].y = bb.min.y;
        poly[1].x = next_x;
        poly[1].y = bb.min.y;
        poly[2].x = next_x;
        poly[2].y = bb.max.y;
        poly[3].x = *x;
        poly[3].y = bb.max.y;
        
        // intersect with this expolygon
        // append results to return value
        polygons_append(*polygons, intersection(poly, to_polygons(*this)));
    }
}

void
ExPolygon::get_trapezoids2(Polygons* polygons, double angle) const {
    ExPolygon clone = *this;
    clone.rotate(PI / 2 - angle, Point(0, 0));
    clone.get_trapezoids2(polygons);
    for (Polygons::iterator polygon = polygons->begin(); polygon != polygons->end(); ++polygon)
        polygon->rotate(-(PI / 2 - angle), Point(0, 0));
}

void
ExPolygon::get_trapezoids3_half(Polygons* polygons, float spacing) const {

    // get all points of this ExPolygon
    Points pp = *this;

    if (pp.empty()) return;

    // build our bounding box
    BoundingBox bb(pp);

    // get all x coordinates
    int min_x = pp[0].x, max_x = pp[0].x;
    std::vector<coord_t> xx;
    for (Points::const_iterator p = pp.begin(); p != pp.end(); ++p) {
        if (min_x > p->x) min_x = p->x;
        if (max_x < p->x) max_x = p->x;
    }
    for (int x = min_x; x < max_x-spacing/2; x += spacing) {
        xx.push_back(x);
    }
    xx.push_back(max_x);
    //std::sort(xx.begin(), xx.end());

    // find trapezoids by looping from first to next-to-last coordinate
    for (std::vector<coord_t>::const_iterator x = xx.begin(); x != xx.end() - 1; ++x) {
        coord_t next_x = *(x + 1);
        if (*x == next_x) continue;

        // build rectangle
        Polygon poly;
        poly.points.resize(4);
        poly[0].x = *x +spacing / 4;
        poly[0].y = bb.min.y;
        poly[1].x = next_x -spacing / 4;
        poly[1].y = bb.min.y;
        poly[2].x = next_x -spacing / 4;
        poly[2].y = bb.max.y;
        poly[3].x = *x +spacing / 4;
        poly[3].y = bb.max.y;

        // intersect with this expolygon
        // append results to return value
        polygons_append(*polygons, intersection(poly, to_polygons(*this)));
    }
}

// While this triangulates successfully, it's NOT a constrained triangulation
// as it will create more vertices on the boundaries than the ones supplied.
void
ExPolygon::triangulate(Polygons* polygons) const
{
    // first make trapezoids
    Polygons trapezoids;
    this->get_trapezoids2(&trapezoids);
    
    // then triangulate each trapezoid
    for (Polygons::iterator polygon = trapezoids.begin(); polygon != trapezoids.end(); ++polygon)
        polygon->triangulate_convex(polygons);
}

void
ExPolygon::triangulate_pp(Polygons* polygons) const
{
    // convert polygons
    std::list<TPPLPoly> input;
    
    ExPolygons expp = union_ex(simplify_polygons(to_polygons(*this), true));
    
    for (ExPolygons::const_iterator ex = expp.begin(); ex != expp.end(); ++ex) {
        // contour
        {
            TPPLPoly p;
            p.Init(int(ex->contour.points.size()));
            //printf(PRINTF_ZU "\n0\n", ex->contour.points.size());
            for (Points::const_iterator point = ex->contour.points.begin(); point != ex->contour.points.end(); ++point) {
                p[ point-ex->contour.points.begin() ].x = point->x;
                p[ point-ex->contour.points.begin() ].y = point->y;
                //printf("%ld %ld\n", point->x, point->y);
            }
            p.SetHole(false);
            input.push_back(p);
        }
    
        // holes
        for (Polygons::const_iterator hole = ex->holes.begin(); hole != ex->holes.end(); ++hole) {
            TPPLPoly p;
            p.Init(hole->points.size());
            //printf(PRINTF_ZU "\n1\n", hole->points.size());
            for (Points::const_iterator point = hole->points.begin(); point != hole->points.end(); ++point) {
                p[ point-hole->points.begin() ].x = point->x;
                p[ point-hole->points.begin() ].y = point->y;
                //printf("%ld %ld\n", point->x, point->y);
            }
            p.SetHole(true);
            input.push_back(p);
        }
    }
    
    // perform triangulation
    std::list<TPPLPoly> output;
    int res = TPPLPartition().Triangulate_MONO(&input, &output);
    if (res != 1) CONFESS("Triangulation failed");
    
    // convert output polygons
    for (std::list<TPPLPoly>::iterator poly = output.begin(); poly != output.end(); ++poly) {
        long num_points = poly->GetNumPoints();
        Polygon p;
        p.points.resize(num_points);
        for (long i = 0; i < num_points; ++i) {
            p.points[i].x = coord_t((*poly)[i].x);
            p.points[i].y = coord_t((*poly)[i].y);
        }
        polygons->push_back(p);
    }
}

void
ExPolygon::triangulate_p2t(Polygons* polygons) const
{
    ExPolygons expp = simplify_polygons_ex(*this, true);
    
    for (ExPolygons::const_iterator ex = expp.begin(); ex != expp.end(); ++ex) {
        // TODO: prevent duplicate points

        // contour
        std::vector<p2t::Point*> ContourPoints;
        for (Points::const_iterator point = ex->contour.points.begin(); point != ex->contour.points.end(); ++point) {
            // We should delete each p2t::Point object
            ContourPoints.push_back(new p2t::Point(point->x, point->y));
        }
        p2t::CDT cdt(ContourPoints);

        // holes
        for (Polygons::const_iterator hole = ex->holes.begin(); hole != ex->holes.end(); ++hole) {
            std::vector<p2t::Point*> points;
            for (Points::const_iterator point = hole->points.begin(); point != hole->points.end(); ++point) {
                // will be destructed in SweepContext::~SweepContext
                points.push_back(new p2t::Point(point->x, point->y));
            }
            cdt.AddHole(points);
        }
        
        // perform triangulation
        cdt.Triangulate();
        std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();
        
        for (std::vector<p2t::Triangle*>::const_iterator triangle = triangles.begin(); triangle != triangles.end(); ++triangle) {
            Polygon p;
            for (int i = 0; i <= 2; ++i) {
                p2t::Point* point = (*triangle)->GetPoint(i);
                p.points.push_back(Point(point->x, point->y));
            }
            polygons->push_back(p);
        }

        for(std::vector<p2t::Point*>::iterator it = ContourPoints.begin(); it != ContourPoints.end(); ++it) {
            delete *it;
        }
    }
}

Lines
ExPolygon::lines() const
{
    Lines lines = this->contour.lines();
    for (Polygons::const_iterator h = this->holes.begin(); h != this->holes.end(); ++h) {
        Lines hole_lines = h->lines();
        lines.insert(lines.end(), hole_lines.begin(), hole_lines.end());
    }
    return lines;
}

std::string
ExPolygon::dump_perl() const
{
    std::ostringstream ret;
    ret << "[" << this->contour.dump_perl();
    for (Polygons::const_iterator h = this->holes.begin(); h != this->holes.end(); ++h)
        ret << "," << h->dump_perl();
    ret << "]";
    return ret.str();
}

BoundingBox get_extents(const ExPolygon &expolygon)
{
    return get_extents(expolygon.contour);
}

BoundingBox get_extents(const ExPolygons &expolygons)
{
    BoundingBox bbox;
    if (! expolygons.empty()) {
        for (size_t i = 0; i < expolygons.size(); ++ i)
			if (! expolygons[i].contour.points.empty())
				bbox.merge(get_extents(expolygons[i]));
    }
    return bbox;
}

BoundingBox get_extents_rotated(const ExPolygon &expolygon, double angle)
{
    return get_extents_rotated(expolygon.contour, angle);
}

BoundingBox get_extents_rotated(const ExPolygons &expolygons, double angle)
{
    BoundingBox bbox;
    if (! expolygons.empty()) {
        bbox = get_extents_rotated(expolygons.front().contour, angle);
        for (size_t i = 1; i < expolygons.size(); ++ i)
            bbox.merge(get_extents_rotated(expolygons[i].contour, angle));
    }
    return bbox;
}

extern std::vector<BoundingBox> get_extents_vector(const ExPolygons &polygons)
{
    std::vector<BoundingBox> out;
    out.reserve(polygons.size());
    for (ExPolygons::const_iterator it = polygons.begin(); it != polygons.end(); ++ it)
        out.push_back(get_extents(*it));
    return out;
}

bool remove_sticks(ExPolygon &poly)
{
    return remove_sticks(poly.contour) || remove_sticks(poly.holes);
}

} // namespace Slic3r
