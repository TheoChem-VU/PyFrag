# kill $(ps aux | grep '[b]okeh' | awk '{print $2}') 2>/dev/null

RESULTDIR="$( pwd -P )"
portnumber=$RESULTDIR/stocks/port.txt
port=$(<$portnumber)


bokeh serve --port $port --show $RESULTDIR/stocks &
