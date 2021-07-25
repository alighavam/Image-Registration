function hdScore = hd(img1,img2,searchWin)

[d12,d21] = surfd(img1,img2,searchWin,0);
hdScore = max([max(d12),max(d21)]);