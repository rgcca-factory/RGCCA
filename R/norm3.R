norm3=function (vec) 
{
  a <- sqrt(sum(vec^2,na.rm=T))
  if (a == 0) 
    a <- 0.05
  return(a)
}