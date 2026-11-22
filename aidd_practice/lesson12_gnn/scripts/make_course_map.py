"""根据系列课程结构绘制第12篇定位图。课程关系示意，不是训练结果图。"""
from common import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
plt.rcParams.update({'font.sans-serif':['Hiragino Sans GB','Arial Unicode MS'],'axes.unicode_minus':False})
def main():
 fig,ax=plt.subplots(figsize=(9,11.5));fig.patch.set_facecolor('#F4F8FA');ax.set_facecolor('#F4F8FA');ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
 def text(x,y,t,fs=12,c='#24495A',weight='normal',ha='left'):
  ax.text(x,y,t,transform=ax.transAxes,fontsize=fs,color=c,fontweight=weight,ha=ha,va='center',linespacing=1.65)
 def box(x,y,w,h,fc='#FFFFFF',ec='#D5E2E5'):
  ax.add_patch(FancyBboxPatch((x,y),w,h,transform=ax.transAxes,boxstyle='round,pad=.012,rounding_size=.012',facecolor=fc,edgecolor=ec,lw=1.1))
 def arrow(y1,y2):ax.add_patch(FancyArrowPatch((.5,y1),(.5,y2),transform=ax.transAxes,arrowstyle='-|>',mutation_scale=12,color='#167F77',lw=1.4))
 text(.03,.971,'「零基础入门 AI 制药」',21,weight='bold')
 text(.03,.932,'第12篇在整套学习计划中的位置',17,c='#0E8075',weight='bold')
 box(.03,.785,.94,.105)
 text(.05,.865,'01 先建立认知、工具与实践底座',14,weight='bold')
 text(.05,.832,'① 起始 → ②③ 分子语言 → ⑤ 工具 → ⑥ 全流程 → ⑧ MVP',11.8)
 text(.05,.803,'④ 政策、⑦ 行业：理解目标与边界，不是技术先修门槛',11,c='#6A7F88')
 arrow(.771,.742)
 text(.03,.719,'02 从跑通 MVP，深入检查各个技术模块',14,c='#0E8075',weight='bold')
 cards=[(.598,'第9篇  结构与口袋','研究对象准备好了吗？结构、可信度与口袋','拿到结构，不等于已经确定结合位点'),(.483,'第10篇  三维空间与对接','分子怎样摆放？姿势与评分分别验证','摆得像，不等于活性能排得准'),(.368,'第11篇  分子表征与活性预测','指纹 + RF：标签、划分与适用范围','随机划分好成绩，不等于新骨架也可靠'),(.253,'第12篇  本篇：图神经网络与表示升级','指纹换成图：增益、波动和代价是什么？','24次拟合，用学习曲线检验是否值得升级')]
 for k,(y,title,body,core) in enumerate(cards):
  current=k==3;box(.05,y,.9,.073,fc='#0D756C' if current else '#FFFFFF',ec='#0D756C' if current else '#D5E2E5')
  text(.07,y+.064,title,14,c='#FFFFFF' if current else '#285164',weight='bold')
  text(.07,y+.037,body,12,c='#F1FFFB' if current else '#355D6A')
  text(.07,y+.012,core,10.4,c='#CBECE5' if current else '#72818A')
  if k<3:arrow(y-.014,y-.029)
 arrow(.238,.211)
 box(.03,.103,.94,.095,fc='#FFF9EE',ec='#E5D6B8')
 text(.05,.174,'03 下一站：第13篇 分子生成（计划）',14,c='#A16D21',weight='bold')
 text(.05,.143,'从给分子打分，到提出新候选；生成高分仍需独立验证',11.5)
 text(.05,.116,'后续规划：ADMET · 逆合成 · 干湿闭环；证据回流升级⑧ MVP',10.5,c='#6F7D83')
 text(.03,.066,'本篇属于“分子性质预测”支线，不是全新独立赛道。',11.5,c='#0E8075',weight='bold')
 text(.03,.046,'计划导航：W7建直觉 · W12—W15练模型 · W17—W18写对照报告（方法迁移）',9.4,c='#0E8075')
 text(.03,.018,'注：箭头表示课程承接；研发中的结构、预测和属性评估可并行或迭代。',9.4,c='#74858D')
 path=ROOT/'figures/图00_第12篇课程定位.png';fig.savefig(path,dpi=170,bbox_inches='tight',facecolor=fig.get_facecolor());plt.close(fig);print(path)
if __name__=='__main__':main()
