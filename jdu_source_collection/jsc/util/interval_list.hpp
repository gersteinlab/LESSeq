#ifndef _jsc_util_interval_hpp_included_
#define _jsc_util_interval_hpp_included_

#include <boost/config.hpp>

#include <math.h>

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>

#include <jsc/util/log.hpp>

using namespace std;
using namespace boost;
using namespace boost::lambda;

namespace jsc
{
namespace util
{

//! A function which checks whether two intervals overlap.
/*!
 * \return whether [a, b) and [c, d) overlap.
 */
template<class T, class Compare>
bool interval_overlap(T const & a, T const & b, T const & c, T const & d, Compare comp = Compare())
{
	return (!(comp(b, c) || comp(d, a)));
}

//! A function which computes the overlapping length of two intervals.
/*!
 * \return the overlapping length of [a, b) and [c, d).
 */
template<class T, class Compare>
double compute_interval_overlap(T const & a, T const & b, T const & c, T const & d, Compare comp = Compare())
{
	double length = 0;
	if (!(comp(b, c) || comp(d, a)))
	{
		if (comp(a, c))
		{
			if (comp(b, d))
			{
				length = b - c;
			}
			else
			{
				length = d - c;
			}
		}
		else
		{
			if (comp(b, d))
			{
				length = b - a;
			}
			else
			{
				length = d - a;
			}
		}
	}
	return length;
}

//! A class representing a list of intervals.
/*!
 * The intervals in the list are sorted according to their start positions and do not overlap with each other.
 */
template <class T, class Compare = less<T> >
class interval_list
{
private:
	//! The comparator.
	Compare comp;
	//! The internal representation of the list of intervals: [starts[0], end[0]), [starts[1], end[1]), ...
	vector<T> starts;
	//! The internal representation of the list of intervals: [starts[0], end[0]), [starts[1], end[1]), ...
	vector<T> ends;

public:
	//! A constructor.
	interval_list()
		: comp(Compare())
	{
	}

	//! A constructor taking an explicit comparator.
	explicit
	interval_list(Compare const & c)
		: comp(c)
	{
	}

	//! Copy constructor.
	interval_list(interval_list<T, Compare> const & il)
		: comp(il.comp), starts(il.starts), ends(il.ends)
	{
	}

	//! Gets the number of intervals.
	unsigned long get_num_intervals() const
	{
		return starts.size();
	}

	//! Gets the starts vector.
	vector<T> const & get_starts() const
	{
		return starts;
	}

	//! Gets the ends vector.
	vector<T> const & get_ends() const
	{
		return ends;
	}

	//! Counts the number of intervals in the list that overlap with the given interval.
	/*!
	 * \sa count_overlap_il().
	 */
	long count_overlap(T const & start, T const & end) const
	{
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return 0;	// always return 0 in this case
		}

		typename vector<T>::const_iterator s_starts_itr, e_ends_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);
		e_ends_itr = lower_bound(ends.begin(), ends.end(), end);

		long s_starts_idx, e_ends_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		e_ends_idx = distance(ends.begin(), e_ends_itr);
		long i;
		long total_count = 0;
		for (i = max((long)0, s_starts_idx - 1); i <= min((long)starts.size() - 1, e_ends_idx); ++i)
		{
			if (interval_overlap<T, Compare>(starts[i], ends[i], start, end, comp))
			{
				total_count++;
			}
		}

		return total_count;
	}

	//! Counts the number of intervals in the list that overlap with the given interval list.
	/*!
	 * \sa count_overlap().
	 */
	long count_overlap_il(interval_list<T, Compare> const & il) const
	{
		long total_count = 0;
		int n = starts.size();

		for (int i = 0; i < n; ++i)
		{
			if (il.check_overlap(starts[i], ends[i]))
			{
				total_count++;
			}
		}
		return total_count;
	}

	//! Checks whether the interval list overlaps with the interval [start, end).
	bool check_overlap(T const & start, T const & end) const
	{
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return false;	// always return false in this case
		}

		typename vector<T>::const_iterator s_starts_itr, e_ends_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);
		e_ends_itr = lower_bound(ends.begin(), ends.end(), end);

		long s_starts_idx, e_ends_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		e_ends_idx = distance(ends.begin(), e_ends_itr);
		long i;
		for (i = max((long)0, s_starts_idx - 1); i <= min((long)starts.size() - 1, e_ends_idx); ++i)
		{
			if (interval_overlap<T, Compare>(starts[i], ends[i], start, end, comp))
			{
				return true;
			}
		}

		return false;
	}

	static interval_list<T> get_interval_overlap(T const & a, T const & b, T const & c, T const & d, Compare comp = Compare())
	{
		interval_list<T> il;
		if (!(comp(b, c) || comp(d, a)))
		{
			if (comp(a, c))
			{
				if (comp(b, d))
				{
					il.add_interval(c, b);
				}
				else
				{
					il.add_interval(c, d);
				}
			}
			else
			{
				if (comp(b, d))
				{
					il.add_interval(a, b);
				}
				else
				{
					il.add_interval(a, d);
				}
			}
		}
		return il;
	}

	//! Gets the overlap of the interval list with the interval [start, end).
	interval_list<T> get_overlap_il(T const & start, T const & end) const
	{
		interval_list<T> il;
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return il;	// always return 0 in this case
		}

		typename vector<T>::const_iterator s_starts_itr, e_ends_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);
		e_ends_itr = lower_bound(ends.begin(), ends.end(), end);

		long s_starts_idx, e_ends_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		e_ends_idx = distance(ends.begin(), e_ends_itr);
		long i;
		for (i = max((long)0, s_starts_idx - 1); i <= min((long)starts.size() - 1, e_ends_idx); ++i)
		{
			il.add_interval_list(
					get_interval_overlap(starts[i], ends[i], start, end, comp));
		}

		return il;
	}

	//! Get the overlap of the interval list with another interval list.
	interval_list<T> get_overlap_il(interval_list<T, Compare> const & il) const
	{
		interval_list<T> rst_il;
		int n = il.starts.size();
		for (int i = 0; i < n; i++)
		{
			rst_il.add_interval_list(get_overlap_il(il.starts[i], il.ends[i]));
		}
		return rst_il;
	}

	//! Computes the overlap of the interval list with the interval [start, end).
	/*!
	 * \sa compute_overlap_il().
	 */
	double compute_overlap(T const & start, T const & end) const
	{
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return 0;	// always return 0 in this case
		}

		typename vector<T>::const_iterator s_starts_itr, e_ends_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);
		e_ends_itr = lower_bound(ends.begin(), ends.end(), end);

		long s_starts_idx, e_ends_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		e_ends_idx = distance(ends.begin(), e_ends_itr);
		long i;
		double total_length = 0;
		for (i = max((long)0, s_starts_idx - 1); i <= min((long)starts.size() - 1, e_ends_idx); ++i)
		{
			total_length += compute_interval_overlap<T, Compare>(starts[i], ends[i], start, end, comp);
		}

		return total_length;
	}

	bool check_single_overlap(T const & start, T const & end, double const & threshold) const
	{
		bool found = false;
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return false;	// always return false in this case
		}

		typename vector<T>::const_iterator s_starts_itr, e_ends_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);
		e_ends_itr = lower_bound(ends.begin(), ends.end(), end);

		long s_starts_idx, e_ends_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		e_ends_idx = distance(ends.begin(), e_ends_itr);
		long i;
		for (i = max((long)0, s_starts_idx - 1); i <= min((long)starts.size() - 1, e_ends_idx); ++i)
		{
			if (threshold <= compute_interval_overlap<T, Compare>(starts[i], ends[i], start, end, comp)) {
				found = true;
				break;
			}
		}

		return found;
	}


	//! Computes the overlap of the interval list with another interval list.
	/*!
	 * \sa compute_overlap().
	 */
	double compute_overlap_il(interval_list<T, Compare> const & il) const
	{
		double total_length = 0;
		int n = il.starts.size();
		for (int i = 0; i < n; i++)
		{
			total_length += compute_overlap(il.starts[i], il.ends[i]);
		}
		return total_length;
	}

	//! Computes the total length of the intervals in this interval list.
	double compute_total_length() const
	{
		double total_length = 0;
		int n = starts.size();
		for (int i = 0; i < n; i++)
		{
			total_length += ends[i] - starts[i];
		}
		return total_length;
	}

	//! Finds in the interval list the index of the interval which contains the given interval [start, end).
	/*!
	 * \return the index of the interval containing [start, end); otherwise starts.size().
	 * \sa contains_interval(T const & start, T const & end)
	 */
	long idx_containing_interval(T const & start, T const & end) const
	{
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return 0;	// always return 0 in this case
		}

		typename vector<T>::const_iterator s_starts_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);

		long s_starts_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		if (s_starts_idx >= 0 && s_starts_idx < (long)starts.size() &&
				!comp(start, starts[s_starts_idx]) &&
				!comp(ends[s_starts_idx], end))
		{
			return s_starts_idx;
		}
		if (s_starts_idx - 1 >= 0 && s_starts_idx - 1 < (long)starts.size() &&
				!comp(start, starts[s_starts_idx - 1]) &&
				!comp(ends[s_starts_idx - 1], end))
		{
			return s_starts_idx - 1;
		}

		return starts.size();
	}

	//! Checks whether the interval list contains [start, end).
	/*!
	 * \sa idx_containing_interval()
	 */
	bool contains_interval(T const & start, T const & end) const
	{
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return true;
		}

		typename vector<T>::const_iterator s_starts_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);

		long s_starts_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		if (s_starts_idx >= 0 && s_starts_idx < (long)starts.size() &&
				!comp(start, starts[s_starts_idx]) &&
				!comp(ends[s_starts_idx], end))
		{
			return true;
		}
		if (s_starts_idx - 1 >= 0 && s_starts_idx - 1 < (long)starts.size() &&
				!comp(start, starts[s_starts_idx - 1]) &&
				!comp(ends[s_starts_idx - 1], end))
		{
			return true;
		}

		return false;
	}

	//! Fill in gaps if their length is less than the given threshold
	void fill_in_gaps(double const & threshold)
	{
		int n = starts.size();
		if (n == 0)
		{
			return;
		}

		vector<T> new_starts, new_ends;
		T last_end;
		new_starts.push_back(starts[0]);
		last_end = ends[0];
		for (int i = 1; i < n; ++i)
		{
			if (starts[i] - last_end < threshold)
			{
				// fill in the gap
			}
			else
			{
				// copy the last end and current start
				new_ends.push_back(last_end);
				new_starts.push_back(starts[i]);
			}
			// update last end
			last_end = ends[i];
		}
		new_ends.push_back(last_end);
		starts = new_starts;
		ends = new_ends;
	}

	//! Add [start, end) to the interval list.
	/*!
	 * Overlapping intervals in the list will be merged.
	 * \sa add_starts_sizes() and add_interval_list()
	 */
	interval_list<T> & add_interval(T const & start, T const & end)
	{
		if (!comp(start, end))	/* [start, end) == nil */
		{
			return *this;
		}

		typename vector<T>::iterator s_starts_itr, s_ends_itr, e_starts_itr, e_ends_itr;
		s_starts_itr = lower_bound(starts.begin(), starts.end(), start);
		s_ends_itr = lower_bound(ends.begin(), ends.end(), start);
		e_starts_itr = lower_bound(starts.begin(), starts.end(), end);
		e_ends_itr = lower_bound(ends.begin(), ends.end(), end);

		int s_starts_idx, s_ends_idx, e_starts_idx, e_ends_idx;
		s_starts_idx = distance(starts.begin(), s_starts_itr);
		s_ends_idx = distance(ends.begin(), s_ends_itr);
		e_starts_idx = distance(starts.begin(), e_starts_itr);
		e_ends_idx = distance(ends.begin(), e_ends_itr);
		bool s_on_interval = (s_starts_idx - s_ends_idx == 1);
		bool e_on_interval = (e_starts_idx - e_ends_idx == 1);

		typename vector<T>::iterator starts_itr = starts.erase(s_starts_itr, e_starts_itr);
		typename vector<T>::iterator ends_itr = ends.erase(s_ends_itr, e_ends_itr);
		if (s_on_interval && e_on_interval)
		{
		}
		else if (!s_on_interval && !e_on_interval)
		{
			starts.insert(starts_itr, start);
			ends.insert(ends_itr, end);
		}
		else if (s_on_interval && !e_on_interval)
		{
			ends.insert(ends_itr, end);
		}
		else
		{
			starts.insert(starts_itr, start);
		}

		return *this;
	}

	//! Add another interval list to the interval list.
	/*!
	 * \sa add_starts_sizes() and add_interval()
	 */
	void add_interval_list(interval_list<T, Compare> const & il)
	{
		int n = il.starts.size();
		for (int i = 0; i < n; i++)
		{
			add_interval(il.starts[i], il.ends[i]);
		}
	}

	//! Add another interval list to the interval list.
	/*!
	 * \param starts the starts vector of the interval list to be added.
	 * \param sizes  the sizes vector of the interval list to be added.
	 * \sa add_interval_list() and add_interval()
	 */
	void add_starts_sizes(vector<T> const & starts, vector<T> const & sizes)
	{
		int n = starts.size();
		for (int i = 0; i < n; i++)
		{
			add_interval(starts[i], starts[i] + sizes[i]);
		}
	}

	//! Add another interval list to the interval list.
	/*!
	 * \param starts the starts vector of the interval list to be added.
	 * \param ends the ends vector of the interval list to be added.
	 * \sa add_interval_list() and add_interval()
	 */
	void add_starts_ends(vector<T> const & starts, vector<T> const & ends)
	{
		int n = starts.size();
		for (int i = 0; i < n; i++)
		{
			add_interval(starts[i], ends[i]);
		}
	}

	//! Checks whether the total coverage of the interval list overlaps with that of another interval list.
	/*!
	 * coverage = [min(starts), (max(ends))
	 */
	bool coverage_overlap(interval_list<T, Compare> const & il2)
	{
		int n1 = starts.size();
		int n2 = il2.starts.size();
		if (n1 == 0 || n2 == 0)
		{
			return false;
		}
		return interval_overlap<T, Compare>(starts[0], ends[n1 - 1], il2.starts[0], il2.ends[n2 - 1], comp);
	}

	T find_min_uncovered(T start_from) const {
		if (get_num_intervals() == 0) {
			return start_from;
		} else {
			typename vector<T>::const_iterator s_starts_itr;
			s_starts_itr = lower_bound(starts.begin(), starts.end(), start_from);

			long s_starts_idx;
			s_starts_idx = distance(starts.begin(), s_starts_itr);
			long i;
			for (i = max((long)0, (long)s_starts_idx - 1);
					i < (long)get_num_intervals(); ++i)
			{
				if (starts[i] > start_from) {
					return start_from;
				} else if (ends[i] > start_from) {
					start_from = ends[i];
				} else {
					continue;
				}
			}
			return start_from;
		}
	}
};

//! Default printing of an interval_list.
template<class T, class Compare>
ostream & operator << (ostream& os, interval_list<T, Compare> const & il)
{
	int n = il.get_starts().size();
	for (int i = 0; i < n; i++)
	{
		os << "[" << il.get_starts()[i] << "," << il.get_ends()[i] << "), ";
	}
	return os;
}


} /* end of util */

} /* end of jsc */

#endif
