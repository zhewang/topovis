-- http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx

SELECT objid,ra,dec,modelmag_r,modelmag_r-5*log10(4.28E+08*z) as abs_mag_r, z
FROM SpecPhoto
WHERE modelmag_r < 17.72
AND z >= 0.03 and z <= 0.09
AND ra >= 120 and ra <= 260
AND dec > 0 and dec <= 60
AND class = "GALAXY"
AND modelmag_r-5*log10(4.28E+08*z) < -19.6 and modelmag_r-5*log10(4.28E+08*z) > -21.1
