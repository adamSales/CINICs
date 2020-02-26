data{
 int N;
 real Y[N];
}
parameters{
 real mu;
}
model{
 for(i in 1:N)
  Y[i]~normal(mu,1);
 mu~normal(0,1);
}

