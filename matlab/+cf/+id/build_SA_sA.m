function [SA, sA] = build_SA_sA()
% Construct SA and sA such that vec(A) = SA*alpha + sA, for A.8 form.
% Dimensions: vec(A) is 64x1, alpha is 22x1.

n = 8;
na = 22;

% Base A0 with ones on diagonal where fixed and zeros elsewhere (excluding alpha slots)
A0 = zeros(n,n);
A0(1,1)=1; A0(2,2)=1; A0(3,3)=1; A0(4,4)=1; A0(5,5)=1; A0(6,6)=1; A0(7,7)=1; A0(8,8)=1;

% SA maps each alpha(k) to its position in A
SA = zeros(n*n, na);

% helper to set mapping for A(r,c) = alpha(k)
    function setpos(r,c,k)
        idx = sub2ind([n n], r, c);
        SA(idx,k) = 1;
    end

% Fill from A.8
setpos(2,1,1);
setpos(3,1,2);
setpos(6,1,3);
setpos(7,1,4);
setpos(8,1,5);

setpos(3,2,6);
setpos(6,2,7);
setpos(7,2,8);
setpos(8,2,9);

setpos(7,3,10);
setpos(8,3,11);

setpos(5,4,12);
setpos(6,4,13);
setpos(7,4,14);
setpos(8,4,15);

setpos(5,6,16);

setpos(6,5,17);
setpos(7,5,18);
setpos(8,5,19);

setpos(7,6,20);
setpos(8,6,21);

setpos(8,7,22);

sA = A0(:);
end
