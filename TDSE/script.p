set grid
set terminal gif animate delay 0.3
unset key
set yrange [0:10]
set output "out.gif"
do for [i=0:5000:100] {
    k = i/10
    set title "t = " . k . "s"
    fname = sprintf("%d.txt", i)
    plot fname using 1:2 with lines
}