function out = mylog2(in)
  out = log2(in);
  out(~in) = 0;

end