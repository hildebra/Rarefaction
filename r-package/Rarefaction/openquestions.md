# Open questions
- is the returned divvs object of the expected size?
- is the retCnts object of the correct size?
- what is the object div? It currently is empty
- what happens in line 78 in IO.cpp '''ss>>ID; ss>>num;''' is 'num' now a float with value 10010... or a float with value cumsum?
- because some samples do not have a to small totSum<dept some samples do not have a retCnt as rarefying is not performed. Currently the retCnt is not linked to the sample count. This must be changed?

- what are we rarefying for? OTUS or species?
