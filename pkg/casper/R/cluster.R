load('casper/data/dm3.newTxs.RData');
trans = lapply(newTxs, as.integer);

for (v in trans)
{
  for (e in v)
  {
    ex2va[e] = 
  }
}

eclust = NULL;
nextid = 1;
while (isEmpty(ex2va) == false)
{
  c = nextid;
  nextid++;
  s = key(ex2va[0]);
  queue = { s };
  while (isEmpty(queue) == false)
  {
    eclust[queue] = c;
    queue = ex2va[queue];
  }
}

return eclust;
