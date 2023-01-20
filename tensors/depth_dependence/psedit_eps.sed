# change just plotted line widths
s/userlinewidth setlinewidth/userlinewidth 3 mul setlinewidth/
# change all line widths
#s/gnulinewidth 5.000/gnulinewidth 15.000/
# change dash length
s/dl {10 mul}/dl {30 mul}/
# change cross size
s/1.000 UP/1.500 UP/
# change bounding box
s/BoundingBox: 50 50/BoundingBox: 50 45/
