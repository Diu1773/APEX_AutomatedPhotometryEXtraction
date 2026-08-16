---
title: 'APEX: Automated Photometry EXtraction for cluster CMDs and multi-night light curves'
tags:
  - Python
  - astronomy
  - photometry
  - light curves
  - color-magnitude diagram
  - isochrone fitting
authors:
  - name: "FIRST LAST"            # TODO: replace with your name
    orcid: 0000-0000-0000-0000    # TODO: add your ORCID
    affiliation: 1
affiliations:
  - name: "Your institution"      # TODO
    index: 1
date: 20 June 2026
bibliography: paper.bib
---

# Summary

`APEX` (Automated Photometry EXtraction) is an open-source Python application for
end-to-end astronomical aperture and PSF photometry. It takes raw FITS frames
through source detection, WCS plate solving, master-catalog construction, and
forced aperture photometry, and then branches into two science modes: cluster
color-magnitude diagrams (CMD) with PARSEC isochrone fitting, and multi-night
differential light curves with detrending and period analysis (Lomb–Scargle,
PDM, BLS). It runs as a PyQt5 desktop application, and as a scriptable
command-line pipeline that needs no GUI toolkit installed: calibration,
detection, astrometry, catalog construction, aperture and PSF photometry,
zero-point calibration and the isochrone fit all live in `apex.analysis`, and
the desktop windows drive those same objects rather than copies of them. The
batch isochrone fit declines to run until the settings that decide the answer —
the colours, the age window, any prior — are written in the configuration, since
library defaults yield a confidently wrong age. The interactive master-ID
editor, the CMD viewer and the light-curve branch remain desktop-only.

`APEX` is built on the Astropy ecosystem [@astropy2022], using `photutils` for
detection and aperture/PSF photometry, `astroquery` for Gaia/SIMBAD access, and
`scipy`/`numpy` for the numerical core. It exports to community submission
formats (AAVSO Extended File Format, ExoClock, ExoFOP/TFOP).

# Statement of need

Amateur and small-observatory astronomers performing variable-star and cluster
photometry rely on tools such as AstroImageJ [@collins2017] and HOPS. These are
strong for single-object transit follow-up but are GUI-centric and offer limited
headless automation, and they do not target cluster CMD/isochrone work or
multi-night ensemble light curves. `APEX` addresses this gap with (1) a unified
desktop + scriptable pipeline in which the window and the batch run are the same
objects, not two implementations kept in agreement, (2) first-class cluster CMD
and isochrone fitting, and (3) multi-night merging and period analysis at scale
(tens to thousands of frames).

# Validation

`APEX`'s photometry is validated *directly* against independent software on real
observations, rather than only through downstream products. Re-measuring its
forced-aperture star positions with SExtractor [@bertin1996] and with IRAF/DAOPHOT
[@stetson1987] at the same aperture agrees with `APEX` to ~3 mmag (robust MAD), with
>99 % of stars within 0.05 mag, and a cross-match against Gaia DR3 synthetic photometry
[@gaiadr3] shows the colours agree to ~10 mmag with no colour-dependent trend. End-to-end
runs on six open and globular clusters (0.85–7.8 kpc) recover published distance moduli to
~4 %. These checks are automated (`apex validate --suite crosscheck`) and documented with
full tables and reproduction commands in the project documentation.

# Acknowledgements

`APEX` builds on Astropy [@astropy2022], photutils, and SExtractor
[@bertin1996].

# References
