#kill $(ps aux | grep '[b]okeh' | awk '{print $2}') 2>/dev/null

RESULTDIR="$( pwd -P )"

bokeh serve  --show $RESULTDIR/stocks &

PID=$!
sleep    $RESULTCHECK/2
kill -9 $PID

exit 0
