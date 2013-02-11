julian_day = function(t)
{
    # Based on http://aa.usno.navy.mil/faq/docs/JD_Formula.php

    return( sapply(t, function(d) 
    {
        d = with_tz(d,"UTC")
        K = year(d)
        M = month(d)
        I = day(d)
        UT = hour(d) + minute(d)/60 + second(d)/60^2

        JD = 367*K - floor( 7*(K+floor((M+9)/12)) / 4) + floor(275 * M / 9) + I + 1721013.5 + UT/24 - 0.5*sign(100*K+M-190002.5) + 0.5

        return(JD)
    }
    ))
}
