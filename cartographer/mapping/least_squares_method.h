#ifndef LEAST_SQUARES_METHOD_H
#define LEAST_SQUARES_METHOD_H

#include <deque>
#include <tuple>
#include <numeric>
#include <cmath>

#include <boost/circular_buffer.hpp>

#include "cartographer/common/time.h"
#include "cartographer/sensor/odometry_data.h"
#include "cartographer/transform/rigid_transform.h"

#include "glog/logging.h"

#define SQR(x)((x)*(x))

namespace cartographer
{

namespace lsm
{

    inline bool isnanvector( const transform::Rigid3d::Vector &v )
    {
        return   std::isnan(v.x())
              && std::isnan(v.y())
              && std::isnan(v.z());
    }

class LeastSquaresMethod
{
public:
    using Vector3D = transform::Rigid3d::Vector;
    using ultra_long_int = long long int;

    template<typename container>
    void initCoefficientsByOdometryData(container data)
    {
        std::tie(_a, _b) = getLSMcoefficients( data );
    }

    void setFirstTime(int64 firstTime)
    {
        _firstTime = firstTime;
    }

    Vector3D getTranslateByTime( const common::Time &time )
    {
        int64 xTime = common::ToUniversal(time) - _firstTime;

        return _a * xTime + _b;
    }

//! Private fields
private:
    Vector3D _a;
    Vector3D _b;
    int64 _firstTime;

//! Private methods
private:
    template<typename container>
    transform::Rigid3d::Vector sumByTranslation(const container &odometryList)
    {
        return std::accumulate( odometryList.begin(), odometryList.end(), transform::Rigid3d::Vector{0.0, 0.0, 0.0},
                                [&]( transform::Rigid3d::Vector left,
                                     const sensor::OdometryData &right ) -> transform::Rigid3d::Vector
                                {
                                    return left + right.pose.translation();
                                });
    }

    template<typename container>
    ultra_long_int sumByTime(const container &odometryList)
    {
        return std::accumulate( odometryList.begin(), odometryList.end(), 0L,
                                [&]( ultra_long_int left,
                                     const sensor::OdometryData &right ) -> ultra_long_int
                                {
                                    return left + (common::ToUniversal(right.time) - _firstTime);
                                });
    }

    template<typename container>
    transform::Rigid3d::Vector sumByTranslationAndTimeMultiplies(const container &odometryList)
    {
        return std::accumulate( odometryList.begin(), odometryList.end(), transform::Rigid3d::Vector{0.0, 0.0, 0.0},
                                [&]( const transform::Rigid3d::Vector &left,
                                     const sensor::OdometryData &right ) -> transform::Rigid3d::Vector
                                {
                                    ultra_long_int time = ToUniversal(right.time) - _firstTime;
                                    return left + time * right.pose.translation();
                                });
    }

    template<typename container>
    ultra_long_int sumByTimeSquare(const container &odometryList)
    {
        return std::accumulate( odometryList.begin(), odometryList.end(), 0L,
                                [&]( ultra_long_int left,
                                     const sensor::OdometryData &right ) -> ultra_long_int
                                {
                                    ultra_long_int time = ToUniversal( right.time ) - _firstTime;
                                    return left + SQR(time);
                                });
    }

    template<typename container>
    std::tuple<Vector3D, Vector3D> getLSMcoefficients(const container &data)
    {
        auto Ex = double( sumByTime( data ) );
        auto Ey = sumByTranslation( data );
        auto Ex_pow_2 = double( sumByTimeSquare( data ) );
        auto Exy = sumByTranslationAndTimeMultiplies( data );

        int64 n = int64( data.size() );
        auto a =  ( n * Exy - Ex * Ey )
                / ( n * Ex_pow_2 - SQR(Ex) );

        auto b = ( Ey - a * Ex ) / n;

        if( isnanvector(a) || isnanvector(b) || true )
        {
            LOG(INFO) << "LSM: a or b nan received:";
            LOG(INFO) << "LSM: a: " << a.x() << " " << a.y() << " " << a.z();
            LOG(INFO) << "LSM: b: " << b.x() << " " << b.y() << " " << b.z();
            LOG(INFO) << "LSM: Coefficients:";
            LOG(INFO) << "LSM: Ex:   " << Ex;
            LOG(INFO) << "LSM: Ey:   " << Ey.x() << " " << Ey.y() << " " << Ey.z();
            LOG(INFO) << "LSM: Ex^2: " << Ex_pow_2;
            LOG(INFO) << "LSM: Exy:  " << Exy.x() << " " << Exy.y() << " " << Exy.z();
            LOG(INFO) << "LSM: n:    " << n;

            LOG(INFO) << "LSM: printing time array:";
            for( auto it = data.begin(); it != data.end(); ++it )
            {
                LOG(INFO) << "LSM: " << common::ToUniversal(it->time) - _firstTime << " | "
                          << it->pose.translation().x() << " "
                          << it->pose.translation().y() << " "
                          << it->pose.translation().z();
            }
        }

        return std::tuple<Vector3D, Vector3D>( a, b );
    }

};

}  // namespace LSF

}  // cartographer


#endif
