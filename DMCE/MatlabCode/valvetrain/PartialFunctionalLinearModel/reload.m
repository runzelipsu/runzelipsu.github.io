clear;
nb=500;
ne=5000;
load 'flm-mse' mse11 mse21 mse12 mse22;
load 'flm-mmse' mmse11 mmse21 mmse12 mmse22;
load 'flm-varmse' varmse11 varmse21 varmse12 varmse22;
load 'flm-ratio' rmse1 rmse2 rvarmse1 rvarmse2;

for i=10:ne
	rmse1(i)=mmse11(i)/mmse21(i);
	rmse2(i)=mmse12(i)/mmse22(i);
	rvarmse1(i)=varmse11(i)/varmse21(i);
	rvarmse2(i)=varmse12(i)/varmse22(i);
	rsmse1(i)=log(mse11(i)/mse21(i))/log(10);
   rsmse2(i)=log(mse12(i)/mse22(i))/log(10);
end;
figure(21)
a=[rsmse1(nb:ne)',rsmse2(nb:ne)'];
BOXPLOT(a);
xlabel('');
ylabel('log_{10}(r_i)');

figure(4)
plot(mmse11(nb:ne));
ylabel('MMSE_0^1');

figure(5)
%plot(mmse(2,1,:));
plot(mmse21(nb:ne));
ylabel('MMSE_0^2');

figure(6)
plot(mmse12(nb:ne));
ylabel('MMSE_1^1');

figure(7)
plot(mmse22(nb:ne));
ylabel('MMSE_1^2');

figure(8)
plot(varmse11(nb:ne));
ylabel('VarMSE_0^1');

figure(9)
plot(varmse21(nb:ne));
ylabel('VarMSE_0^2');

figure(10)
plot(varmse12(nb:ne));
ylabel('VarMSE_1^1');

figure(11)
plot(varmse22(nb:ne));
ylabel('VarMSE_1^2');

figure(12)
plot(mse11(nb:ne));
ylabel('MSE_0^1');

figure(13)
%plot(mmse(2,1,:));
plot(mse21(nb:ne));
ylabel('MSE_0^2');

figure(14)
plot(mse12(nb:ne));
ylabel('MSE_1^1');

figure(15)
plot(mse22(nb:ne));
ylabel('MSE_1^2');

figure(16)
plot(rmse1(nb:ne));
ylabel('RMSE^1');

figure(17)
plot(rmse2(nb:ne));
ylabel('RMSE^2');

figure(18)
plot(rvarmse1(nb:ne));
ylabel('RVARMSE^1');

figure(19)
plot(rvarmse2(nb:ne));
ylabel('RVARMSE^2');



figure(20)
subplot(2,2,1);
plot([nb:1:ne],rmse1(nb:ne));
%title('a');
axis([nb ne  10.5 11.5]);
xlabel('k');
ylabel('MMSE^{(1)}_0/MMSE^{(2)}_0');
subplot(2,2,2);
plot([nb:1:ne],rmse2(nb:ne));
axis([nb ne  14 16]);
xlabel('k');
ylabel('MMSE^{(1)}_1/MMSE^{(2)}_1');
subplot(2,2,3);
plot([nb:1:ne],rvarmse1(nb:ne));
axis([nb ne  8 11]);
xlabel('k');
ylabel('VarMSE^{(1)}_0/VarMSE^{(2)}_0');
subplot(2,2,4);
plot([nb:1:ne],rvarmse2(nb:ne));
axis([nb ne  22 35]);
xlabel('k');
ylabel('VarMSE^{(1)}_1/VarMSE^{(2)}_1');
print -dpsc flm-ratio.eps;

