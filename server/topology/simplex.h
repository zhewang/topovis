#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>

#include "metric_space.h"
#include "../geometry/point_incs.h"


class Simplex  {
	public:
		Simplex();
		Simplex(const std::vector<int> & _simplex, MetricSpace* _metricSpace);
		Simplex(const std::vector<int> & _simplex, Points &points, std::map<int,int> loc_map);
    Simplex(const std::vector<int> & _simplex, double _distance, std::map<int,int> loc_map);
		~Simplex();

		int dim() const { return simplex.size()-1; }
		int vertex(int _i) const { return simplex[_i]; }
    int location(int _i) const { return loc_map.at(_i); }
		int min_vertex() const { return *std::min_element(simplex.begin(), simplex.end()); }
		double get_simplex_distance() const { return cached_distance; }
    std::vector<int> as_vector() { return simplex;};

    //MetricSpace* get_metric_space() { return this->metric_space; }

		std::string id() const {
			return uid;
		}

		std::vector<Simplex> faces(double** distances=NULL);

		inline bool operator==(const Simplex& _simp) const {
            return this->simplex == _simp.simplex;
        }

		inline bool operator<(const Simplex& _simp) const {
			int our_dim = this->dim(), other_dim = _simp.dim();
			// lower dimensional simplices come before higher dimensional ones
			/*
			if(our_dim < other_dim)
				return true;
			else if(other_dim < our_dim)
				return false;
				*/

			// if dimensions are 0, then just use index for comparison
			if(our_dim==0 && other_dim==0)
				return simplex[0] < _simp.vertex(0);

			// otherwise, take largest edge distance on simplex for comparison
			double our_edge_dist = this->get_simplex_distance(), other_edge_dist = _simp.get_simplex_distance();
            if(our_edge_dist < other_edge_dist)
                return true;
            if(other_edge_dist < our_edge_dist)
                return false;

			// if they are equivalent, then lower dimensional simplices precede higher dimensional ones
            if(our_dim < other_dim)
                return true;
            if(our_dim > other_dim)
                return false;

			// resolve equal distance, equal dim
            for(int i = 0; i < other_dim+1; i ++) {
                if(this->vertex(i) < _simp.vertex(i))
                    return true;
                if(this->vertex(i) > _simp.vertex(i))
                    return false;
            }
            return false;
		}

		friend std::ostream& operator <<(std::ostream &out, const Simplex & _simplex)  {
			for(int i = 0; i <= _simplex.dim(); i++)  {
                if( i > 0) { out << " "; }
				out <<  _simplex.vertex(i) << "(" << _simplex.location(_simplex.vertex(i)) << ")";
			}
			return out;
		}

		void compute_simplex_distance(double ** _distances)  {
			for(unsigned i = 0; i < simplex.size(); i++)  {
				for(unsigned j = 0; j < i; j++)  {
					double next_dist = _distances[loc_map[simplex[i]]][loc_map[simplex[j]]];
					cached_distance = next_dist > cached_distance ? next_dist : cached_distance;
				}
			}
		}

	private:
		void compute_simplex_distance()  {
			for(unsigned i = 0; i < simplex.size(); i++)  {
				for(unsigned j = 0; j < i; j++)  {
					double next_dist = metric_space->distance(loc_map[simplex[i]], loc_map[simplex[j]]);
					cached_distance = next_dist > cached_distance ? next_dist : cached_distance;
				}
			}
		}

		void compute_simplex_distance_from_points(Points &points) {
      for(unsigned i = 0; i < simplex.size(); i++)  {
        for(unsigned j = 0; j < i; j++)  {
          double next_dist = (points[loc_map[simplex[i]]] - points[loc_map[simplex[j]]]).l2Norm();
          cached_distance = next_dist > cached_distance ? next_dist : cached_distance;
        }
      }
    }

    void compute_uid() {
      char unique_id[10*(simplex.size()+1)];
      sprintf(unique_id, "%u", simplex[0]);
      for(unsigned i = 1; i < simplex.size(); i++) {
        sprintf(unique_id, "%s-%u", unique_id, simplex[i]);
      }
      uid = std::string(unique_id);
    }

    // <PointID, PositionID>, after subdivision, some inserted points should be
    // treated as the same PointID, but they have different positions.
    std::map<int, int> loc_map;

		std::vector<int> simplex;
		MetricSpace* metric_space;
		double cached_distance;
    std::string uid;
};


struct global_compare {
    static std::map<std::string, int> order_map;

    bool operator() (const Simplex& left, const Simplex& right) const{
        return global_compare::order_map.at(left.id()) < global_compare::order_map.at(right.id());
    }
};

struct lex_compare {
    bool operator() (const Simplex& left, const Simplex& right) const{
        int left_dim = left.dim(), right_dim = right.dim();

        //if dimensions are 0, then just use index for comparison
        if(left_dim==0 && right_dim==0)
            return left.vertex(0) < right.vertex(0);

        //lower dimensional simplices precede higher dimensional ones
        if(left_dim < right_dim)
            return true;
        if(left_dim > right_dim)
            return false;

        for(int i = 0; i < left_dim+1; i ++) {
            if(left.vertex(i) < right.vertex(i))
                return true;
            if(left.vertex(i) > right.vertex(i))
                return false;
        }
        return false;
    }
};

#endif
