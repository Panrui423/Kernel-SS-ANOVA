load("ADNIdatatrain.Rdata")
load("ADNIdatatest.Rdata")

Y.train = log(data.train$Score + 2) # logrithm transformation 
N = length(Y.train)
xfun0.train = data.train$Xfun
s.train = data.train$Xscalar  # c(Intercept,  gender, Edu, age);  normalized age, Edu


Y.test = log(data.test$Score + 2)
xfun0.test = data.test$Xfun
s.test = data.test$Xscalar


Y.all = c(c(Y.train),c(Y.test))
s = rbind(s.train, s.test)
Xfun <- rbind(xfun0.train, xfun0.test)

m = length(Y.all)

label_left = read.csv('label-left-12.csv')
label_right = read.csv('label-right-12.csv')
