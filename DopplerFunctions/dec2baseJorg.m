%produce a vector of numbers from index in base 10

function outvec = dec2baseJorg(index, newbase)
indcopy = index;
finished = false;
outvec = zeros( 1, floor( log(index)/log(newbase) )+1 ); 
ctr = 1;
while ~finished,
    outvec(ctr) = mod( indcopy, newbase);
    indcopy = (indcopy-outvec(ctr) )/newbase;
    if indcopy == 0,
        finished = true;
    end
    ctr = ctr + 1;
end
outvec = fliplr( outvec);