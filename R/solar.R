library(lubridate)
library(timezone)

refraction_correction = function(elev)
{
    if (is.na(elev))
        return(NA)

    elev_rad = elev * pi / 180

    if (elev > 85) {
        r = 0
    } else if (elev > 5) {
        r = 58.1/tan(elev_rad) - 0.07/tan(elev_rad)^3 + 0.000086/tan(elev_rad)^5
    } else if (elev > -0.575) {
        r = 1735 + elev*(-518.2+elev*(103.4+elev*(-12.79+elev*0.711)))
    } else {
        r = -20.772/tan(elev_rad)
    }

    return(r / 3600)
}

solar = function(t, x, y, p4s="", tz="")
{
    stopifnot(is.POSIXt(t))
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    stopifnot(is.character(p4s))
    stopifnot(is.character(tz))

    #stopifnot(all(lapply(list(t, x, y, p4s, tz), length) == 1))
    stopifnot(all(lapply(list(x, y, p4s, tz), length) == 1))

    if (tz == "")
        tz = find_tz(x,y,p4s=p4s)

    #attr(t,"tzone") = tz
    t = force_tz(t,tz)
    t = with_tz(t,"UTC")

    jd = julian_day(t)
    jc = (jd-2451545) / 36525    

    
    K = year(t)
    M = month(t)
    I = day(t)
    UT = hour(t) + minute(t)/60 + second(t)/60^2


    crds = cbind(x,y)
    if (p4s != "")
        crds = project(crds, p4s, inverse=TRUE)


    lat = crds[,2]
    lat_rad = lat * pi / 180
    long = crds[,1]
    long_rad = long * pi / 180

    geom_mean_long_sun = (280.46646 + jc * (36000.76983 + jc * 0.0003032)) %% 360
    geom_mean_anom_sun = 357.52911 + jc * (35999.05029 - 0.0001537 * jc)
    eccent_earth_orbit = 0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)

    sun_eq_of_ctr = sin((pi/180) * geom_mean_anom_sun) * (1.914602-jc*(0.004817+0.000014*jc)) +
                    sin((pi/180) * 2 * geom_mean_anom_sun) * (0.019993-0.000101*jc)           +
                    sin((pi/180) * 3 * geom_mean_anom_sun) * 0.000289

    sun_true_long = geom_mean_long_sun + sun_eq_of_ctr
    sun_true_anom = geom_mean_anom_sun + sun_eq_of_ctr

    sun_rad_vec = (1.000001018 * (1-eccent_earth_orbit^2)) / (1 + eccent_earth_orbit*cos((pi/180)*sun_true_anom))
    sun_app_long = sun_true_long - 0.00569 - 0.00478*sin( (pi/180)*(125.04 - 1934.136*jc) )

    mean_obliq_ecliptic = 23 + (26 + ((21.448 - jc*(46.815 + jc*(0.00059 - jc*0.001813))))/60)/60
    obliq_corr = mean_obliq_ecliptic + 0.00256 * cos( (pi/180) * (125.04 - 1934.136*jc) )

    sun_rt_ascen = (180/pi) * atan2( cos((pi/180) * obliq_corr) * sin((pi/180) * sun_app_long), cos((pi/180) * sun_app_long) )
    sun_declin = (180/pi) * asin( sin((pi/180) * obliq_corr) * sin((pi/180) * sun_app_long) )
    
    var_y = tan((pi/180) * obliq_corr/2)^2

    eq_of_time = 4 * (180/pi) * (var_y*sin(2*(pi/180)*geom_mean_long_sun) - 2*eccent_earth_orbit*sin((pi/180)*geom_mean_anom_sun) + 4*eccent_earth_orbit*var_y*sin((pi/180)*geom_mean_anom_sun)*cos(2*(pi/180)*geom_mean_long_sun) - 0.5*var_y^2*sin(4*(pi/180)*geom_mean_long_sun)-1.25*eccent_earth_orbit^2*sin(2*(pi/180)*geom_mean_anom_sun))

    true_solar_time = (UT/24 * 1440 + eq_of_time + 4*long) %% 1440

    hour_angle = true_solar_time/4 - sign(true_solar_time)*180


    solar_zenith = (180/pi) * acos(sin(lat_rad) * sin( (pi/180) * sun_declin ) + cos(lat_rad) * cos( (pi/180) * sun_declin ) * cos( (pi/180) * hour_angle ))
    solar_elev = 90 - solar_zenith
    
    approx_atmo_refrac = sapply(solar_elev, refraction_correction)

    solar_elev_corr = solar_elev + approx_atmo_refrac


    tmp = (180/pi) * (acos( ( (sin(lat_rad) * cos(pi/180 * solar_zenith)) - sin(pi/180 * sun_declin) ) / ( cos(lat_rad) * sin(pi/180 * solar_zenith))))
    solar_azimuth = ifelse(hour_angle > 0, 180 + tmp, 540 - tmp) %% 360

    ha_sunrise = (180/pi) * acos( cos(90.833 * pi/180) / (cos(lat_rad) * cos(sun_declin * pi/180))
                                 -tan(lat_rad) * tan(sun_declin * pi/180)
                                )

    solar_noon = t
    second(solar_noon) = 0
    minute(solar_noon) = 0
    hour(solar_noon) = 12

    solar_noon = force_tz(solar_noon - seconds(round(60*(4*long+eq_of_time))), tz)

    sunrise = solar_noon - seconds(round(60*4*ha_sunrise))
    sunset  = solar_noon + seconds(round(60*4*ha_sunrise))


    return(
        data.frame(
            declin = sun_declin,
            zenith = solar_zenith,
            elev = solar_elev,
            elev_corr = solar_elev_corr,
            azimuth = solar_azimuth,
            noon = solar_noon,
            sunrise = sunrise,
            sunset = sunset
        )
    )  
}
