
function res = fb_gamlike(P,K)

a=P(1);
b=P(2);

res = -sum( log( fb_gampdf(K, a, b) ) );

 	  	 

