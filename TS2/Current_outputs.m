function [ Ik1, Ito, Ikr, Iks, ICaL, INaK, INa, INaCa, If, IbNa, IbCa, Istim, ICaT, IpCa, Irel,Iup,Ileak ] = Current_outputs( Time, values, model_parameter_inputs )
Ik1 = zeros(size(Time));
Ito = zeros(size(Time));
Ikr = zeros(size(Time));
Iks = zeros(size(Time));
ICaL = zeros(size(Time));
INaK = zeros(size(Time));
INa = zeros(size(Time));
INaCa = zeros(size(Time));
IpCa = zeros(size(Time));
If = zeros(size(Time));
IbNa = zeros(size(Time));
IbCa = zeros(size(Time));
Irel = zeros(size(Time));
Iup = zeros(size(Time));
Ileak = zeros(size(Time));
Istim = zeros(size(Time));
ICaT = zeros(size(Time));



for i= 1:size(values,1)
    [temp data] =  ipsc_function(Time(i), values(i,:), model_parameter_inputs);
    Ik1(i) = data(1);
    Ito(i) = data(2);
    Ikr(i) = data(3);
    Iks(i) = data(4);
    ICaL(i) = data(5);
    INaK(i) = data(6);
    INa(i) = data(7);
    INaCa(i) = data(8);
    IpCa(i) = data(9);
    If(i) = data(10);
    IbNa(i) = data(11);
    IbCa(i) = data(12);
    Irel(i) = data(13);
    Iup(i) = data(14);
    Ileak(i) = data(15);
    Istim(i) = data(16);
    ICaT(i)= data(17);
   
end


end

