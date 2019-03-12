
function res = fb_gamfit(R)

avg = mean(R);

a=nmsmax( @gamfit_search, 1, [], [], avg, R );

b=a/avg;      

res=[a 1/b];


function res = gamfit_search( a, avg, R )

b=a/avg;      

res = -fb_gamlike([a 1/b], R);

 	  	 

