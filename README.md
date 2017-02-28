# network-model
1.BA_SIR.m 
  是借用大牛的http://modelinginfectiousdiseases.org/ 对其进行了修改。
  可实现BA，随机网络，小世界，spatial网络的SIR和SIS传播

2.BA_SIR_commuinity(单词拼错==)
  在1的基础上加入社区，使用matlab自带的Kmeans函数，划分10个社区，每个社区人数不同，可为0，可为1
  实现  if 结点n：S->I
          n所在的社区的S的Rate+1，I-0.01（即如果社区中有人感染了，那么还没感染的人接触率就会变大，已经感染的人的恢复率就会减小）
        if 结点n：I->R
          n所在的社区的S的Rate-1，I+0.01（即如果社区中有人恢复了，那么还没感染的人接触率就会减小，已经感染的人的恢复率就会增大）

3.BA_SIR_com_noK
  在1的基础上加入社区，同等划分10个社区，其余与2相同

4.BA_SIR_com_diffS
  S的状态分为主动和被动，主动：设置为总结点数的1/10，改变接触率，暂时设为10，可调整
