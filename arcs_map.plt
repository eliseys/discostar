set terminal qt	 

set view equal xyz

set cbrange [0.0:18.0]

splot "wire.data" u 1:2:3 w l lw 0.1 lc "black" notitle,\
      "arcs.data" u 1:2:3:4 palette pt 7 ps 0.5 notitle,\
      "arcs.data" u ($4 == 15 ? $1 : 1/0):2:3:4 pt 7 ps 1.0 lt rgb "navy" notitle,\
      "arcs.data" u ($4 == 16 ? $1 : 1/0):2:3:4 pt 7 ps 1.0 lt rgb "red" notitle,\
      "arcs.data" u ($4 == 17 ? $1 : 1/0):2:3:4 pt 7 ps 1.0 lt rgb "gray" notitle,\
      "arcs.data" u ($4 == 18 ? $1 : 1/0):2:3:4 w l lw 0.2 lc "red" notitle

pause -1
#quit
