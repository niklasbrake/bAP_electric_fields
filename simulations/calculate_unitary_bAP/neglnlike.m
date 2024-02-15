function output = neglnlike(B1,X0,X1,y)
    model = X0 + B1*X1;
    output = sum(log(model) + y./model);
end