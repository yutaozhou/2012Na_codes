def cal_logg(vmag, a_v, plx, teff, mass, feh=0.0, test=False):
    mbol_sun, logg_sun, teff_sun = 4.74, 4.44, 5777.
    x = np.log10(teff) - 3.52
    if (np.log10(teff) >= 3.50) & (np.log10(teff) <= 3.67) & (feh >= -0.50) & (feh <= 0.20) or \
       (np.log10(teff) >= 3.56) & (np.log10(teff) <= 3.67) & (feh >= -1.50) & (feh <= -0.50) or \
       (np.log10(teff) >= 3.58) & (np.log10(teff) <= 3.67) & (feh >= -2.50) & (feh <= -1.50) or \
       (np.log10(teff) >= 3.61) & (np.log10(teff) <= 3.67) & (feh >= -3.00) & (feh <= -2.50):
        bc = -0.05531 / x - 0.6177 + 4.420 * x - 2.669 * x**2 + 0.6943 * x * feh - 0.1071 * feh - 0.008612 * feh**2
    elif (np.log10(teff) >= 3.65) & (np.log10(teff) <= 3.96) & (feh >= -0.50) & (feh <= 0.20) or \
         (np.log10(teff) >= 3.65) & (np.log10(teff) <= 3.83) & (feh >= -1.50) & (feh <= -0.50) or \
         (np.log10(teff) >= 3.65) & (np.log10(teff) <= 3.80) & (feh >= -2.50) & (feh <= -1.50) or \
         (np.log10(teff) >= 3.65) & (np.log10(teff) <= 3.74) & (feh >= -3.00) & (feh <= -2.50):
        bc = -0.09930 / x + 0.02887 + 2.275 * x - 4.425 * x**2 + 0.3505 * x * feh - 0.05558 * feh - 0.005375 * feh**2
    else:
        raise ValueError('out of ranges of applications!')
    if test:  # 2014Bergemann_A&A_565_89
        print 'The test correction in A_V should be for SFD98'
        e_bv = a_v / 3.1
        if e_bv > 0.1:
            e_bv = 0.035 + 0.65*e_bv
            a_v = 3.1*e_bv
            mbol_star = vmag + bc + 5 * np.log10(plx) + 5.0 - a_v  # plx in arcsec
    else:
        mbol_star = vmag + bc + 5 * np.log10(plx) + 5.0 - a_v  # plx in arcsec
    logL = -0.4*(mbol_star - mbol_sun)
    logg = logg_sun + np.log10(mass) + 4 * np.log10(teff / teff_sun) + 0.4 * (mbol_star - mbol_sun)
    return '{:.3f}'.format(logg), '{:.3f}'.format(logL)