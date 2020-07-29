# Orbits-Project-Results
The following is the result of the pipeline I have built to compute the Ephemeris of an object in the Solar System. The pipeline takes three sets observations(each  set containing a time of observation, right ascension and declination), begining date, time interval and number of outputs required. 

The finalpipeline.py contains the complete pipeline.
The finalpipeline.ipynb is the Jupyter NB of the above pipeline. You can easily access the notebook from NBViewer.

Here is the result of an attempt to compute the Orbital Elements and Ephemeris from real data from JPL's HORIZONS.

HORIZONS object - [Asteroid (2019 PH3)](https://ssd.jpl.nasa.gov/horizons.cgi#results)


OBITAL ELEMENTS:

-----------------

latus rectum = 1.6443376922928608

eccentricity = 0.6134240790325921  (99.79% accurate)

a = 2.636377998582081  (99.64% accurate)

Period = 1563.541110708  (99.46% accurate)

ArgP = 173.8753979327009 (99.99% accurate)

Node = 143.78120214112022 (99.98% accurate)

inclination = 15.724593583226332 (99.91% accurate)

TP = 2019-08-12T05:18:02.126 (99.82% acurate)

MA = 0.8703193343309066 (99.99% accurate)


EPHEMERIS:


|         Date        |      RA        |       DEC       |     Distance(AU)    |
|---------------------|----------------|-----------------|---------------------|
| 2019-09-05 00:00:00 | 03h18m09.8912s | -33d54m36.2111s | 0.13246958856451996 |
| 2019-09-25 00:00:00 | 03h10m25.0618s |  -33d09m23.881s |  0.2676040874322413 |
| 2019-10-15 00:00:00 | 02h47m35.4727s | -30d59m17.6802s | 0.40868196411984414 |
| 2019-11-04 00:00:00 | 02h27m24.647s  | -26d42m09.2122s |  0.577543725148535  |
| 2019-11-24 00:00:00 | 02h19m12.0027s | -20d47m24.8923s |  0.7928037612422921 |
| 2019-12-14 00:00:00 | 02h23m12.935s  | -14d19m43.4238s |  1.0610788424173423 |
| 2020-01-03 00:00:00 | 02h36m22.0528s | -08d09m41.8175s |  1.3777573453305743 |
| 2020-01-23 00:00:00 | 02h55m44.851s  |  -02d41m36.929s |  1.7314931971321073 |
| 2020-02-12 00:00:00 | 03h19m17.1919s | +01d56m51.8215s |  2.1077597465570226 |
| 2020-03-03 00:00:00 | 03h45m35.0122s | +05d45m29.6483s |  2.4907088549820346 |
