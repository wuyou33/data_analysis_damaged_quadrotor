function y = butterworth(x,order,wc)

[B,A] = butter(order,wc);
y1 = filtfilt(B,A,flip(x));
y = filtfilt(B,A,flip(y1));
end