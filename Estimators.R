#=========Empirical Study==========

#===pop 1 (situation 1)====

N=5000
N1=4500
N2=500

X=rnorm(N,5,10)
Y=X+rnorm(N,0,1)
x=X+rnorm(N,1,3)
y=Y+rnorm(N,1,3)

X2=rnorm(N2,5,10)
Y2=X2+rnorm(N2,0,1)
x2=X2+rnorm(N2,1,3)
y2=Y2+rnorm(N2,1,3)


Xbar=mean(X)
Ybar=mean(Y)
varX=var(X)
varY=var(Y)


X2bar=mean(X2)
Y2bar=mean(Y2)
varX2=var(X2)
varY2=var(Y2)

V=x-X
U=y-Y

V2=x2-X2
U2=y2-Y2


Vbar=mean(V)
Ubar=mean(U)
varV=var(V)
varU=var(U)
V2bar=mean(V2)
U2bar=mean(U2)
varV2=var(V2)      
varU2=var(U2)

corrXY=cor(X,Y)
corrX2Y2= cor(X2,Y2)


n=550


W2= N2/N
h=2 #4,8
theta= (W2*(h-1))/n
lambda2= (1/n)-(1/N)


#======= Hansen & Hurwitz======


varhh= ((lambda2)*(varY))+((theta)*(varY2))+((lambda2)*(varU))+((theta)*(varU2))


varhhh= ((lambda2)*(varY))+((theta)*(varY2))


#========CV Notations=========

c2Y=(varY)/((Ybar)^2)
c2X=(varX)/((Xbar)^2)
c2Y2=(varY2)/((Ybar)^2)  #(varY2)/((Y2bar)^2)
c2X2=(varX2)/((Xbar)^2) #(varX2)/((X2bar)^2)
c2U=(varU)/((Ubar)^2)
c2V=(varV)/((Vbar)^2)
c2U2=(varU2)/((Ubar)^2)  #(varU2)/((U2bar)^2)
c2V2=(varV2)/((Vbar)^2)   #(varV2)/((V2bar)^2)

sdY= sd(Y)
sdX= sd(X)
sdU= sd(U)
sdV= sd(V)
sdY2= sd(Y2)
sdX2= sd(X2)
sdU2= sd(U2)
sdV2 = sd(V2)

cY= (sdY)/(Ybar)
cX= (sdX)/(Xbar)
cX2= (sdX2)/(Xbar)  #cX2= (sdX2)/(X2bar)
cY2= (sdY2)/(Ybar)    #cY2= (sdY2)/(Y2bar)
cXcY= cX*cY
cY2cX2=cY2*cX2

Nh= ((N+n)/(N-n))
Nhh= ((N+n)/(N-n))^2

Nh4= (1/4)*((N+n)/(N-n))
Nhh4= (1/4)*((N+n)/(N-n))^2

N2= ((N+(2*n))/(N-n))
N22= ((N+(2*n))/(N-n))^2

N24= (1/4)*((N+(2*n))/(N-n))
N224= (1/4)*((N+(2*n))/(N-n))^2


s2U= varU/(Ybar)^2
s2V= varV/(Xbar)^2
s2U2= varU2/(Ybar)^2
s2V2= varV2/(Xbar)^2

#===========Existing Estimators==============

#=====Cochran (1977) t1======


mset1= ((lambda2*(Ybar)^2)*(c2Y + c2X - (2*corrXY*cY*cX))) +

((theta*(Ybar)^2)*(c2Y2+ c2X2- (2*corrX2Y2*cY2*cX2))) +

((lambda2*(Ybar)^2)*(s2U + s2V)) + 

((theta*(Ybar)^2)*(s2U2 + s2V2))

#====t1 without ME====


mset11= ((lambda2*(Ybar)^2)*(c2Y + c2X - (2*corrXY*cY*cX))) +

((theta*(Ybar)^2)*(c2Y2+ c2X2- (2*corrX2Y2*cY2*cX2)))


#======Rao's (1986) t2=========

mset2= ((lambda2*(Ybar)^2)*(c2Y + c2X - (2*corrXY*cY*cX))) +

(theta*varY2) +

((lambda2*(Ybar)^2)*(s2U + s2V)) + 

(theta*varU2)

#====t2 without ME====

mset22= ((lambda2*(Ybar)^2)*(c2Y + c2X - (2*corrXY*cY*cX))) +

(theta*varY2)


#========Singh & Kumar (2008) t3========

mset3= ((lambda2*(Ybar)^2)*(c2Y + c2X - (4*corrXY*cY*cX))) +

((theta*(Ybar)^2)*(c2Y2+ c2X2- (2*corrX2Y2*cY2*cX2))) +

((lambda2*(Ybar)^2)*(s2U + (4*s2V))) + 

((theta*(Ybar)^2)*(s2U2 + s2V2))


#====t3 without ME====

mset33= ((lambda2*(Ybar)^2)*(c2Y + c2X - (4*corrXY*cY*cX))) +

((theta*(Ybar)^2)*(c2Y2+ c2X2- (2*corrX2Y2*cY2*cX2))) 




#========Suggested Estimators=========

#====Azeem Hanif tn1=====

msetn1=((lambda2*(Ybar)^2)*(c2Y+(Nhh*c2X)-(2*Nh*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(Nhh*c2X2)-(2*Nh*corrX2Y2*cY2*cX2))) + 
((lambda2*(Ybar)^2)*(s2U+(Nhh*s2V))) +
((theta*(Ybar)^2)*(s2U2+(Nhh*s2V2)))

#====tn1 without ME====

msetn11=((lambda2*(Ybar)^2)*(c2Y+(Nhh*c2X)-(2*Nh*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(Nhh*c2X2)-(2*Nh*corrX2Y2*cY2*cX2)))


#======tn2=======


msetn2=((lambda2*(Ybar)^2)*(c2Y+(N224*c2X)-(N2*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(N224*c2X2)-(N2*corrX2Y2*cY2*cX2))) + 
((lambda2*(Ybar)^2)*(s2U+(N224*s2V))) +
((theta*(Ybar)^2)*(s2U2+(N224*s2V2)))

#====tn2 without ME====

msetn22=((lambda2*(Ybar)^2)*(c2Y+(N224*c2X)-(N2*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(N224*c2X2)-(N2*corrX2Y2*cY2*cX2))) 



#===optimum value of mu===


mu= (2*((lambda2*corrXY*cY*cX)+(theta*corrX2Y2*cY2cX2)))/((lambda2*((varX+varV)/(Xbar)^2))+(theta*((varX2+varV2)/(Xbar)^2)))

mu2= (mu)^2

mu24= (1/4)*((mu)^2)

#========tn3========


msetn3=((lambda2*(Ybar)^2)*(c2Y+(mu24*c2X)-(mu*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(mu24*c2X2)-(mu*corrX2Y2*cY2*cX2))) + 
((lambda2*(Ybar)^2)*(s2U+(mu24*s2V))) +
((theta*(Ybar)^2)*(s2U2+(mu24*s2V2)))

#====tn3 without ME====

msetn33=((lambda2*(Ybar)^2)*(c2Y+(mu24*c2X)-(mu*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(mu24*c2X2)-(mu*corrX2Y2*cY2*cX2))) 


#============MINE==============


mseha= ((lambda2*(Ybar)^2)*(c2Y+(Nhh4*c2X)-(Nh*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(Nhh4*c2X2)-(Nh*corrX2Y2*cY2*cX2))) + 
((lambda2*(Ybar)^2)*(s2U+(Nhh4*s2V))) +
((theta*(Ybar)^2)*(s2U2+(Nhh4*s2V2))) 

#====Mine without ME====

msehaa= ((lambda2*(Ybar)^2)*(c2Y+(Nhh4*c2X)-(Nh*corrXY*cY*cX))) +
((theta*(Ybar)^2)*(c2Y2+(Nhh4*c2X2)-(Nh*corrX2Y2*cY2*cX2)))



#========Efficiency Comparison==========


ehh= (varhh/varhh)*100
ehhh= (varhhh/varhhh)*100

et1= (varhh/mset1)*100
et11= (varhhh/mset11)*100

et2= (varhh/mset2)*100
et22= (varhhh/mset22)*100

et3= (varhh/mset3)*100
et33= (varhhh/mset33)*100

etn1 = (varhh/msetn1)*100
etn11 = (varhhh/msetn11)*100

etn2 = (varhh/msetn2)*100
etn22=(varhhh/msetn22)*100

etn3 = (varhh/msetn3)*100
etn33 = (varhhh/msetn33)*100

eha = (varhh/mseha)*100
ehaa = (varhhh/msehaa)*100

