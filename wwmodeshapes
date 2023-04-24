function[Modeshapes]=WWmodeshapes(w,nodesdiv,sup1,sup2,Ldiv,E,G,k,I,Ip,A,Rho,L,Mu,c11,c22,c33)
DFGstiff=zeros(3*nodesdiv,3*nodesdiv); %Global dynamic stiffness matrix of the cracked beam.
DFstiff=zeros(3*nodesdiv,3*nodesdiv); 
crackstiff=[1/c11 0 0 -1/c11 0 0; 0 1/c22 0 0 -1/c22 0; 0 0 1/c33 0 0 -1/c33; -1/c11 0 0 1/c11 0 0; 0 -1/c22 0 0 1/c22 0; 0 0 -1/c33 0 0 1/c33]; %Dynamic stiffness matrix of the crack element
for i=1:9 %9 beam elements on the left side of the crack
    Alpha2(i)=Rho*A*w^2*Ldiv(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*Ldiv(i)^4/(E*I);
	rr(i)=I/(A*Ldiv(i)^2);
	ss(i)=E*I/(k*A*G*Ldiv(i)^2);
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	Lambda(i)=(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	Eta(i)=Z(i)/(Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-cos(Phi(i))*cosh(Lambda(i)))+(1-Eta(i)^2)*sin(Phi(i))*sinh(Lambda(i)));
	a1(i)=E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	d1(i)=E*I/Ldiv(i)^3*bb(i)*Tau(i)*(cos(Phi(i))*sinh(Lambda(i))+Eta(i)*sin(Phi(i))*cosh(Lambda(i)))/(Lambda(i)*Phi(i));
	d2(i)=E*I/Ldiv(i)^2*Z(i)*Tau(i)*((Phi(i)+Eta(i)*Lambda(i))*sin(Phi(i))*sinh(Lambda(i))-(Lambda(i)-Eta(i)*Phi(i))*(1-cos(Phi(i))*cosh(Lambda(i))))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=E*I/Ldiv(i)*Tau(i)*(sin(Phi(i))*cosh(Lambda(i))-Eta(i)*cos(Phi(i))*sinh(Lambda(i)));
	d4(i)=-E*I/Ldiv(i)^3*bb(i)*Tau(i)*(sinh(Lambda(i))+Eta(i)*sin(Phi(i)))/(Lambda(i)*Phi(i));
	d5(i)=E*I/Ldiv(i)^2*Z(i)*Tau(i)*(cosh(Lambda(i))-cos(Phi(i)));
	d6(i)=E*I/Ldiv(i)*Tau(i)*(Eta(i)*sinh(Lambda(i))-sin(Phi(i)));
	DFstiff((3*i-2):(3*(i+1)),(3*i-2):(3*(i+1)))=DFstiff((3*i-2):(3*(i+1)),(3*i-2):(3*(i+1)))+[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
end
for i=10:22 %13 beam elements on the right side of the crack
    Alpha2(i)=Rho*A*w^2*Ldiv(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*Ldiv(i)^4/(E*I);
	rr(i)=I/(A*Ldiv(i)^2);
	ss(i)=E*I/(k*A*G*Ldiv(i)^2);
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	Lambda(i)=(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	Eta(i)=Z(i)/(Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-cos(Phi(i))*cosh(Lambda(i)))+(1-Eta(i)^2)*sin(Phi(i))*sinh(Lambda(i)));
	a1(i)=E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	d1(i)=E*I/Ldiv(i)^3*bb(i)*Tau(i)*(cos(Phi(i))*sinh(Lambda(i))+Eta(i)*sin(Phi(i))*cosh(Lambda(i)))/(Lambda(i)*Phi(i));
	d2(i)=E*I/Ldiv(i)^2*Z(i)*Tau(i)*((Phi(i)+Eta(i)*Lambda(i))*sin(Phi(i))*sinh(Lambda(i))-(Lambda(i)-Eta(i)*Phi(i))*(1-cos(Phi(i))*cosh(Lambda(i))))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=E*I/Ldiv(i)*Tau(i)*(sin(Phi(i))*cosh(Lambda(i))-Eta(i)*cos(Phi(i))*sinh(Lambda(i)));
	d4(i)=-E*I/Ldiv(i)^3*bb(i)*Tau(i)*(sinh(Lambda(i))+Eta(i)*sin(Phi(i)))/(Lambda(i)*Phi(i));
	d5(i)=E*I/Ldiv(i)^2*Z(i)*Tau(i)*(cosh(Lambda(i))-cos(Phi(i)));
	d6(i)=E*I/Ldiv(i)*Tau(i)*(Eta(i)*sinh(Lambda(i))-sin(Phi(i)));
	DFstiff((3*i+1):(3*(i+2)),(3*i+1):(3*(i+2)))=DFstiff((3*i+1):(3*(i+2)),(3*i+1):(3*(i+2)))+[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
end
DFGstiff=DFstiff;
DFGstiff(28:33,28:33)=DFstiff(28:33,28:33)+crackstiff;
%Assembling the global dynamic stiffness matrix of the cracked beam.
if isempty(sup1)&&isempty(sup2) %In the case of no supports or restraints for the whole frame
    dis(1:(3*nodesdiv-length(sup1)-length(sup2)-1))=(inv(DFGstiff(1:(3*nodesdiv-length(sup1)-length(sup2)-1),1:(3*nodesdiv-length(sup1)-length(sup2)-1))))*DFGstiff(1:(3*nodesdiv-length(sup1)-length(sup2)-1),(3*nodesdiv-length(sup1)-length(sup2)));
	dis(3*nodesdiv-length(sup1)-length(sup2))=-1; %The highest numbered degree of freedom is assumed to have a unit displacement (or rotation)
	Modeshapes=dis;
else
    DFGstiff([71],:)=[];
	DFGstiff(:,[71])=[];
	DFGstiff([sup1],:)=[];
	DFGstiff(:,[sup1])=[];
	dis(1:(3*nodesdiv-length(sup1)-length(sup2)-1))=(inv(DFGstiff(1:(3*nodesdiv-length(sup1)-length(sup2)-1),1:(3*nodesdiv-length(sup1)-length(sup2)-1))))*DFGstiff(1:(3*nodesdiv-length(sup1)-length(sup2)-1),(3*nodesdiv-length(sup1)-length(sup2)));
	dis(3*nodesdiv-length(sup1)-length(sup2))=-1;
	Modeshapes=dis;
end
end