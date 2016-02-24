#include <iostream>

#include <set>

#include <boost/optional.hpp>
#include <boost/assert.hpp>
#include <boost/log/trivial.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Segment_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck ;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick ;


template < typename K >
class Length_indexed_polyline_2 {
public:
    typedef CGAL::Point_2< K > Point_2 ;
    typedef CGAL::Vector_2< K > Vector_2 ;
    typedef CGAL::Segment_2< K > Segment_2 ;

    typedef std::vector< double >::const_iterator abscisse_iterator ;


    template < typename Iterator >
    Length_indexed_polyline_2( Iterator begin, Iterator end ):
        _points(begin,end)
    {
        BOOST_ASSERT( ! _points.empty() );

        _abscisses.reserve(_points.size());
        _abscisses.push_back(0.0);
        for ( size_t i = 1; i < _points.size(); i++ ){
            const double & last = _abscisses.back();

            Segment_2 segment = Segment_2( _points[i-1], _points[i] );
            _abscisses.push_back(
                last +
                CGAL::sqrt( CGAL::to_double( segment.squared_length() ) )
            ) ;
        }

        BOOST_ASSERT(_abscisses.size() == _points.size() );
    }

    double length() const {
        if ( _abscisses.empty() ){
            return 0.0 ;
        }
        return _abscisses.back();
    }

    Point_2 interpolate( const double & abscisse ) const {
        int index = find_segment(abscisse) ;
        Segment_2 segment = _segment( index ) ;

        Vector_2 st = segment.target() - segment.source() ;

        double k = ( abscisse - _abscisses[index] ) / CGAL::sqrt( CGAL::to_double( segment.squared_length() ) );
        return segment.source() + st * k ;
    }


    abscisse_iterator begin_abscisses() const {
        return _abscisses.begin() ;
    }
    abscisse_iterator end_abscisses() const {
        return _abscisses.end() ;
    }

private:
    std::vector< Point_2  > _points ;
    std::vector< double >   _abscisses ;



    Segment_2 _segment(int index) const {
        return Segment_2(
            _points[index],
            _points[index+1]
        ) ;
    }

    /*
     * find segment index for abscisse
     * @return the segment index -1 if not found
     */
    int find_segment( const double & abscisse ) const {
        if ( abscisse < 0 ){
            return -1 ;
        }
        for ( size_t i = 1; i < _abscisses.size(); i++ ){
            if ( abscisse >= _abscisses[i-1] && abscisse <= _abscisses[i] ){
                return static_cast<int>(abscisse);
            }
        }
        return -1 ;
    }
} ;


template < typename K >
class Polyline_morphing_2 {
public:
    typedef CGAL::Point_2< K >           Point ;
    typedef CGAL::Segment_2< K >         Segment ;
    typedef std::vector< Point >         Polyline ;
    typedef Length_indexed_polyline_2<K> Length_indexed_polyline ;

    template < typename Iterator >
    Polyline_morphing_2(
        Iterator beginSource, Iterator endSource,
        Iterator beginTarget, Iterator endTarget
    ):
        _source(beginSource,endSource),
        _target(beginTarget,endTarget)
    {
        insert_normalized_abscisses(_source);
        insert_normalized_abscisses(_target);
    }

    std::vector< Segment > build_transform_segments() const {
        double sourceLength = _source.length();
        double targetLength = _target.length();

        std::vector< Segment > result ;
        for (
            std::set< double >::const_iterator it = _normalizedAbscisses.begin();
            it != _normalizedAbscisses.end(); ++it
        ){
            result.push_back(Segment(
                _source.interpolate( (*it) * sourceLength ),
                _target.interpolate( (*it) * targetLength )
            ));
        }
        return result;
    }

private:
    Length_indexed_polyline_2<K> _source ;
    Length_indexed_polyline_2<K> _target ;
    std::set< double > _normalizedAbscisses ;


    void insert_normalized_abscisses( const Length_indexed_polyline & lil ) {
        double length = lil.length();
        for (
            typename Length_indexed_polyline::abscisse_iterator it = lil.begin_abscisses();
            it != lil.end_abscisses(); ++it
        ){
            _normalizedAbscisses.insert( *it / length ) ;
        }
    }

} ;


int main( int argc, char *argv[] ){
    std::vector< Epick::Point_2  > source ;
    source.push_back(Epick::Point_2(0.0,0.0));
    source.push_back(Epick::Point_2(1.0,0.0));
    source.push_back(Epick::Point_2(1.0,1.0));

    std::vector< Epick::Point_2  > target ;
    target.push_back(Epick::Point_2(0.0,5.0));
    target.push_back(Epick::Point_2(0.0,6.0));

    Length_indexed_polyline_2< Epick > lil(
        source.begin(),
        source.end()
    ) ;
    std::cout << lil.length() << std::endl ;

    Polyline_morphing_2< Epick > morphing(
        source.begin(), source.end(),
        target.begin(), target.end()
    );
    std::vector< Epick::Segment_2 > segments = morphing.build_transform_segments();

    typedef std::vector< Epick::Segment_2 >::const_iterator segment_iterator ;
    for ( segment_iterator it = segments.begin(); it != segments.end(); ++it ){
        const Epick::Segment_2 & segment = (*it) ;
        std::cout << segment.source().x() << " " << segment.source().y() ;
        std::cout << " -> " ;
        std::cout << segment.target().x() << " " << segment.target().y() ;
        std::cout << std::endl;
    }

//  the max length provides a majorant for hausdorff and frechet distance
//  0 0 -> 0 5
//  1 0 -> 0 5.5
//  1 1 -> 0 6

    return 0;
}
