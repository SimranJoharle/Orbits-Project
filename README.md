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

|         Date        |      RA        |       DEC       |    Distance(AU)    |
|---------------------|----------------|-----------------|--------------------|
| 2020-07-27 00:00:00 | 07h13m18.8216s | +14d14m04.0877s | 4.2179013028567125 |
| 2020-07-28 00:00:00 |  07h14m33.56s  | +14d12m12.7367s | 4.218477917520705  |
| 2020-07-29 00:00:00 | 07h15m47.9696s | +14d10m18.6735s |  4.21887369016094  |
| 2020-07-30 00:00:00 | 07h17m02.047s  | +14d08m21.9427s | 4.219089016943061  |
| 2020-07-31 00:00:00 | 07h18m15.7891s | +14d06m22.5885s | 4.219124247492298  |
| 2020-08-01 00:00:00 | 07h19m29.1929s | +14d04m20.6543s | 4.218979672589334  |
| 2020-08-02 00:00:00 | 07h20m42.2553s |  +14d02m16.183s | 4.218655518551143  |
| 2020-08-03 00:00:00 | 07h21m54.9734s | +14d00m09.2168s | 4.2181519483353105 |
| 2020-08-04 00:00:00 | 07h23m07.3439s | +13d57m59.7975s | 4.217469068720077  |
| 2020-08-05 00:00:00 | 07h24m19.3636s | +13d55m47.9668s | 4.216606942126371  |

