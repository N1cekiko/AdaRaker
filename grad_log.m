function g=grad_log(y,f)


lt=exp(-y*f);

g=-y*lt/(1+lt);