function asdScore = asd(img1,img2,searchWin)

[d12,d21] = surfd(img1,img2,searchWin,0);
asdScore = mean([d12;d21]);