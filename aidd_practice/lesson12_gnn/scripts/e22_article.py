"""只消费已完成实验的JSON/CSV生成文章，不在HTML里编造实验数字。"""
from common import *
import html,re,base64
from urllib.parse import urlsplit

def code(s):return '<pre class="code" style="background:#111111;color:#F1F1F1;padding:16px;border-radius:5px;font-family:Menlo,Consolas,monospace;font-size:12px;line-height:1.75;white-space:pre-wrap;overflow-wrap:anywhere;">'+html.escape(s.strip())+'</pre>'
def fig(n,name,caption):return f'<figure><img src="figures/{name}" alt="图{n}"><figcaption>图 {n}｜{caption}</figcaption></figure>'
def supfig(letter,name,caption):return f'<figure><img src="figures/{name}" alt="补图{letter}"><figcaption>补图 {letter}｜{caption}</figcaption></figure>'
def table(headers,rows):return '<div class="table-wrap"><table><thead><tr>'+''.join('<th>'+str(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in r)+'</tr>' for r in rows)+'</tbody></table></div>'
def main():
 summary=json.loads((ROOT/'results/summary.json').read_text());audit=json.loads((ROOT/'results/e18_audit.json').read_text());b=pd.read_csv(ROOT/'results/benchmark.csv');assert len(b)==24
 full=audit['n_fit'];stats=summary['statistics'];lookup={(r['n_train'],r['model']):r for r in stats};rf=lookup[full,'RF'];gnn=lookup[full,'GNN'];delta=gnn['MAE_mean']-rf['MAE_mean'];winner='GNN' if delta<0 else '指纹＋RF';other='指纹＋RF' if delta<0 else 'GNN';better='低于' if delta<0 else '高于';sizes=sorted(b.n_train.unique());gaps=[lookup[n,'GNN']['MAE_mean']-lookup[n,'RF']['MAE_mean'] for n in sizes];cross=any(a*c<0 for a,c in zip(gaps,gaps[1:]));curve_statement='采样点之间出现了均值排名变化，但不能把两点间的连线交点当作精确门槛。' if cross else '在这四个采样规模上，没有出现两模型均值排名的反转；不能据此断言更大数据下也不会交叉。'
 # 训练成本对资源分配条件敏感，明确只是实测墙钟。
 ratio=gnn['fit_seconds_mean']/rf['fit_seconds_mean'];sens=json.loads((ROOT/'results/e20_sensitivity.json').read_text());full_gnn=b[(b.n_train==full)&(b.model=='GNN')].set_index('seed')
 attachment=json.loads((ROOT/'results/e23_attachment_audit.json').read_text());extra=attachment['scaffold_summary'];paired=attachment['paired_cluster_bootstrap'];correction=attachment['scaling_correction'];erlo=next(r for r in correction['molecules'] if r['tag']=='erlotinib')
 readout=json.loads((ROOT/'results/e25_case_readout.json').read_text())
 headline='把2048位指纹换成一张图，GNN就能赢吗？'
 styles='''body{margin:0;background:#F1F3F5;color:#3F4549;font-family:-apple-system,BlinkMacSystemFont,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif;}main{max-width:680px;margin:auto;background:white;padding:28px 22px 50px;box-sizing:border-box;font-size:15px;line-height:1.85;}h1{font-size:27px;line-height:1.45;color:#153D34;margin:12px 0 18px;}h2{font-size:20px;color:#0F755B;border-left:4px solid #0F755B;padding-left:10px;margin:34px 0 15px;}h3{font-size:16px;color:#2D668E;margin-top:24px;}p{margin:15px 0;}strong{color:#1C3330;}a{color:#08765B;}blockquote,.note{background:#F0F7F4;border-left:3px solid #16866A;margin:20px 0;padding:14px 16px;}blockquote p{margin:0;}figure{margin:25px 0;}figure img{display:block;width:100%;height:auto;}figcaption{font-size:12px;color:#777;text-align:center;line-height:1.65;margin-top:9px;}.code{background:#111;color:#F1F1F1;padding:16px;border-radius:5px;font-family:Menlo,Consolas,monospace;font-size:12px;line-height:1.75;white-space:pre-wrap;overflow-wrap:anywhere;}code{background:#F0F3F5;color:#A34B36;padding:2px 4px;border-radius:3px;}.table-wrap{overflow-x:auto;margin:18px 0;}table{width:100%;border-collapse:collapse;font-size:12px;line-height:1.65;}th,td{border:1px solid #DDE3E2;padding:9px 7px;text-align:left;vertical-align:top;word-break:break-word;}th{background:#EDF5F1;color:#173F35;}.eyebrow{font-size:13px;color:#0F755B;}.meta{font-size:12px;color:#777;}.toc{background:#F5F9F7;border:1px solid #DBE8E1;padding:16px 18px;margin:22px 0;border-radius:5px;}.toc a{text-decoration:none;display:block;margin:6px 0;font-size:14px;}.directory{background:#F4F7F8;color:#45545A;padding:14px;font-size:11px;line-height:1.8;white-space:pre-wrap;overflow-wrap:anywhere;}@media(max-width:450px){main{padding:20px 16px;}h1{font-size:24px;}table{font-size:11px;}th,td{padding:7px 5px;}}'''
 sections=[('s1','一、从固定指纹到可学习图表示：升级了什么？'),('s2','二、一条键怎样“传话”：D-MPNN的最小直觉'),('s3','三、先把比赛条件对齐：10468不等于训练10468'),('s4','四、跑一个真实基线：标准化、训练与验证'),('s5','五、核心成绩单：谁更好，稳定吗？'),('s6','六、学习曲线与成本：交叉点真的出现了吗？'),('s7','七、打开黑盒：GNN到底在看哪里？'),('s8','八、什么时候该用GNN：把三笔账接回决策'),('reproduce','配套目录、复跑命令与资料来源'),('series-index','系列文章索引')]
 out=[f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AI制药第12篇｜{headline}</title><style>{styles}</style></head><body><main>',f'<p class="eyebrow">零基础入门 AI 制药 · 第12篇｜图神经网络与分子表征升级</p><h1>{headline}</h1><p class="meta">从分子指纹到图神经网络 · 附实践代码与配图 · 预计阅读35—40分钟</p>']
 out+=['''<p class="eyebrow">给第一次读这个系列的朋友：先用30秒认识问题</p><p><strong>我们想预测一个分子的实验活性，好决定哪些候选值得优先验证。</strong>把分子画成一张图交给AI，会比先转成一串0和1，预测得更准吗？</p>
<p>上一篇，我们走的是后一条路：先把分子结构转成2048个0和1，这串数字叫作<strong>分子指纹</strong>；再让随机森林从已知实验数据中学习，寻找结构与活性之间的关系。</p>
<p>这条路能用，却也有取舍。指纹保留了许多局部结构模式，但没有完整保留它们之间的连接关系。于是，一个想法很自然：<strong>能不能把原子和化学键组成的图交给模型，让它自己学习哪些结构信息更有用？</strong>这就是本篇要尝试的图神经网络，简称GNN。</p>
<p>看起来，图里保留的连接信息更完整，模型又能自己学特征，应该更准才对。但换一种表示方式，也换了一种学习方法：<strong>多拿到的信息能不能用好，仍然要靠实验检验。</strong></p>
<p>所以这次，我们让指纹模型和图模型学习同一批分子，再用同一批未参与训练的分子检查预测；随后逐步增加训练数据，并重复运行，观察误差、结果波动和训练成本怎样变化。</p>
<blockquote><p>先说一个实测结果：<strong>在本次固定设置下，指纹＋随机森林的平均预测误差更小。</strong>把分子表示换成图，没有自动换来更准的预测。是否值得升级，要把精度、运行波动和成本放在一起衡量。</p></blockquote>
<p>下面先把本篇放回整个学习计划，再从一张分子图开始，看看模型究竟怎样“传话”，这场比较又该怎样做。</p>''']
 out+=['<nav class="toc"><strong>文章结构目录</strong><section id="article-reasoning" style="margin:12px 0 16px;padding:12px 14px;background:#EAF4EF;border-radius:4px;"><p style="margin:0 0 8px;font-weight:700;color:#155D4B;">文章思路</p><p style="margin:0;font-size:14px;line-height:1.9;">先看<strong>为什么要把指纹换成图</strong>，再弄清<strong>图模型怎样学习</strong>；接着把比较规则对齐，用实跑成绩和学习曲线检查<strong>是否更准、需要多少代价</strong>，最后把结论接回MVP与后续分子生成。</p><p style="margin:10px 0 0;font-size:13px;line-height:1.9;color:#47695E;">提出升级问题 → 看懂分子图与消息传递 → 对齐比较规则 → 跑出成绩 → 检查学习曲线与成本 → 审视解释边界 → 决定是否升级</p></section><a href="#course-position">导读：本篇在整个学习计划中的位置</a>'+''.join(f'<a href="#{i}">{t}</a>' for i,t in sections)+'<p style="margin:12px 0 0;font-size:13px;color:#647571;">补充实验入口</p><a href="#attachment-protocol-audit">附件协议复核与配对差值区间</a><a href="#scaling-audit">同一权重的缩放还原纠错</a></nav>']
 out+=['''<h2 id="course-position">导读：本篇在整个学习计划中的位置</h2>
<p>把系列拉远一点看，我们正在经历两个阶段：<strong>先建立全局地图、跑通一个最小流程，再拆开其中的技术模块，逐个检查它们的输入、输出和可信边界。</strong>第12篇属于第二阶段，落在“分子性质预测”支线，是第11篇的继续深入。</p>
<p><strong>①—⑧搭底座：</strong>从为什么学习AI制药、读懂分子语言，到准备工具、认识全流程，最后在第8篇亲手跑通端到端MVP。④政策篇与⑦行业篇帮助理解研发目标和产业边界，不是技术学习的严格先修门槛。</p>
<p><strong>⑨—⑫拆模块：</strong>第9篇检查蛋白结构与口袋，第10篇检查姿势预测与对接评分；第11篇转向“从分子结构预测活性”，建立指纹＋RF基线和评价规则。<strong>本篇接着问：在同一套评价框架里，换成可学习的图表示，是否值得升级？</strong>它没有另开一条赛道，而是在检查MVP中“活性预测”这一站的可选实现。</p>
<figure><img src="figures/图00_第12篇课程定位.png" alt="第12篇在AI制药课程计划中的位置"><figcaption>导读图｜参照系列课程结构更新到第12篇：本篇位于分子性质预测支线；第13篇分子生成及其余主题均按规划标注。箭头表示课程承接，不代表研发流程必须串行执行。</figcaption></figure>
<p>阅读时可以带着三个问题：<strong>承接第11篇，我们保留哪些评估规则？回到第8篇，怎样决定替换哪个预测模块？走向分子生成，怎样防止生成器只是在迎合评分模型？</strong>若尚未跑过MVP，可以先理解本篇概念和结果，再回到配套代码实践；不必把所有支线都当成进入下一阶段的硬门槛。</p>''']
 out+=['''<h3>跟着“每周4小时”计划读，可以怎样对照？</h3>
<p>系列文章按问题展开，学习计划按能力逐步积累，篇号不等于周次。已对照本地《AI制药学习计划_每周4小时版》核对以下导航：</p>''',table(['计划位置','本篇对应内容','怎样使用'],[
['W7：GNN直觉','分子图、消息传递与读出','先读第二节与邻域图，再回到计划中的通用GNN示例'],
['W12—W13：MPNN精读与GNN实战','D-MPNN数据流、特征形状与代码','本篇提供Chemprop路线；不能替代手写MPNN和原计划的ESOL/Tox21任务'],
['W14—W15：划分对照与Chemprop回归','沿用第11篇的评估规则、核对目标缩放','先核对分子索引和量纲，再讨论模型表现'],
['W17—W18：对比报告与双基线','多种子、学习曲线、误差和成本报告','迁移对照方法。计划W18是TDC任务上的ECFP+GBDT，本篇是EGFR上的RF，两者不是同一实验'],
['W29—W34：分子生成','判别模型能回答什么、不能回答什么','将适用范围和误差检查带入后续候选生成与评分']
]),'''<p>更实用的读法是：概念不熟先补W7，代码卡住回到W12—W15，想形成作品则参照W17—W18写报告。无需把本文当作整套计划的“答案册”。</p>''']
 def h(i):out.append(f'<h2 id="{sections[i][0]}">{sections[i][1]}</h2>')
 h(0)
 out+=['''<p>先看第11篇的主角Morgan指纹：围绕每个原子，按照预先定义的规则提取不同半径的局部环境，再哈希、折叠成固定长度向量。<strong>规则是固定的，随机森林学习的是这些模式与活性之间的统计关系。</strong>指纹不是完全没有全局信息——全分子的局部环境会一起进入模型——但连接关系、出现次数与三维构象并没有被本篇位指纹完整保留。</p>
<p>图模型换了一个入口。原子作为节点，键作为边；节点携带原子类型、形式电荷、局部连接及手性等编码，边携带键类型等信息。模型的训练目标不仅影响最后的回归器，也影响中间怎样聚合局部信息。</p>''',table(['表示','提取方式','本篇需要警惕的边界'],[['Morgan位指纹','固定局部环境规则＋哈希折叠','不同环境可能落到同一位；当前配置不编码手性或片段次数'],['2D分子图＋D-MPNN','节点/键初始特征＋可学习消息传递','完整输入连接不等于无损利用全局结构；传播深度和读出会限制信息'],['两者共同缺少什么','本次均未输入三维坐标或蛋白结构','不直接回答结合姿势、构象分布或真实作用机制']]),'''<p>还有一个常被省略的细节：<strong>输入里有手性标记，与模型正确学到手性相关活性差异，是两件事。</strong>本次用布洛芬R/S做特征检查：不编码手性的Morgan相似度为1.000，编码手性后为0.750；当前Chemprop特征器为这对分子生成的节点矩阵不同。它证明“入口能区分”，不能证明“预测一定正确”。</p>''']
 h(1)
 out+=['''<p>可以把消息传递想成：一个局部环境先拿到自身特征，再从相邻环境接收信息，更新自己的向量；重复若干步后，把这些局部表示汇总成一个分子向量，最后预测pIC50。</p>
<p>但我们用的不是任意“通用GNN”，而是Chemprop的<strong>有向键消息传递网络D-MPNN</strong>。它把消息放在方向键上；更新一条方向键时，聚合来源原子周围其他方向的消息，排除立即从反方向回传的那条边。随后将方向键信息汇总到原子，再通过读出得到分子表示。不要把它画成一个实际上不存在的“注意力网络”。<a href="https://pubmed.ncbi.nlm.nih.gov/31361484/">方法来源：Yang等，2019</a>。</p>''',fig(1,'图01_分子图与邻域.png','厄洛替尼同一个NH氮原子的拓扑邻域：0、1、2、3跳分别覆盖1、3、7、12个重原子。橙色固定标中心，绿色标当前邻域内其余原子；各幅分子朝向一致。此图不是D-MPNN的逐层计算图。'),'''<p><strong>这里的“看见”，只说中心原子的局部拓扑范围。</strong>例如3跳的12/29，并不表示模型只读取了12个原子：分子级读出仍会汇总所有原子的表示。图中的跳数按最短键路径计算，不应直接解释成Chemprop的depth参数或真实注意力权重。</p>
<p>这次厄洛替尼图中有29个重原子（即非氢原子）、31条化学键；NH标记中的氢没有单独展开为节点。每条化学键拆成两个方向，所以图中有62条方向边；节点特征矩阵为29×72，方向边特征矩阵为62×14。<strong>方向边数翻倍，不意味着多了一倍化学键。</strong></p>''',code('''# scripts/e18_prepare.py：真实Chemprop接口
from rdkit import Chem
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
mol = Chem.MolFromSmiles(smiles)
graph = SimpleMoleculeMolGraphFeaturizer()(mol)
print(graph.V.shape)  # 厄洛替尼：(29, 72)
print(graph.E.shape)  # 两个方向：(62, 14)'''),'''<p>层数控制信息传播计算的深度，却不是“远处信息会无限流过来”的保证。初始特征、方向键更新、读出和训练目标一起决定模型能利用什么；没有显式输入坐标，这仍然是2D模型。</p>''']
 h(2)
 out+=['''<h3>先过四项审计，再读成绩单</h3>
<p>“公平”不是一句口号，而是一组可以回查的记录。下面把本次已核对的项目，与仍然存在的不对称放在一起：</p>''',table(['审计项','本次如何核对','结论与边界'],[
['数据与划分','固定输入哈希；逐次保存训练索引；检查拟合/验证/测试骨架互不重叠','已核对。两模型同一规模、同一种子使用相同训练分子和1770个测试分子'],
['重复与选择','种子42/43/44全部报告；GNN保存每轮验证轨迹并恢复最优验证权重','已记录。仅一种外层骨架划分，不代表更换测试骨架后的波动'],
['评价口径','MAE/RMSE/R²/Spearman共用实现；从保存的逐分子预测复算指标','已核对。pIC50尺度一致，CSV往返的浮点差异在容差内'],
['配置与预算','正式实验固定预设、未搜索超参数；记录模型配置、轮次和计时口径','已披露，并非完全等预算：RF固定500树，GNN用验证集选轮次，线程数和计时范围也不同']
]),'''<p>这四项不是“盖章即公平”的充分条件。模型输入、学习方式和选择规则仍有差别；它们的作用，是让读者知道<strong>这次比较支持什么，也不支持什么</strong>。</p>''']
 out+=['''<p>第一次公平比较，不从调参开始，而是从分子索引开始。第11篇骨架划分的外层训练集是8698个分子，测试集是1770个；这次从8698个中，再按骨架抽出1550个验证分子，剩下7148个作实际拟合。</p>
<p><strong>验证集负责选轮次，测试集负责汇报表现。</strong>所以学习曲线的最后一个点只能是7148，不能写成“10468个分子全部训练”。后者会把验证和测试都算进去。</p>''',fig(2,'图02_对照协议.png','本次固定数据协议。每个种子只打乱拟合池生成一条嵌套序列；同一规模的RF与GNN读取完全相同的分子索引。'),table(['比较项目','本次设置'],[['模型输入','同一批标准化SMILES、同一份聚合pIC50；RF使用固定ECFP4位矩阵'],['外层划分','直接复用第11篇骨架索引；未换测试集'],['内部验证','GroupShuffleSplit按骨架留出20%组，seed42；产生7148拟合／1550验证'],['学习曲线','500／2000／5000／7148个拟合分子；每种子嵌套随机子集'],['随机种子','42、43、44；固定同一验证和测试集'],['RF配置','500棵树，max_features=1.0，叶节点最少1样本'],['D-MPNN配置','depth3，隐藏维128，2层FFN，dropout0.05，batch128'],['验证与预算','无超参数搜索；GNN最多50轮，至少10轮，8轮无提升停止并恢复最佳验证MAE权重']]),'''<p><strong>这里的“随机种子”，可以理解为复现随机过程的编号。</strong>换一个编号，就能用另一组初始化、数据顺序或抽样重新运行。三个编号都跑，是为了观察结果是否依赖某一次随机设置；它们并不代表重新换了三份测试集。</p>
<p>这里刻意不用“完全等预算”这个说法：RF固定500棵树，GNN通过验证集选择轮次，训练资源和计算路径本来就不同。<strong>我们做的是同数据、固定预设配置的基准对照，不是证明两模型都已调到最优。</strong></p>
<p>当前环境只有CPU。通路检查先跑了300隐藏维模型，之后在正式实验开始前固定128维配置，以完成多种子学习曲线；300维通路结果单独存放，不混进正式表格。因此，本篇GNN与旧版408601参数模型也不是同一个配置。</p>''']
 h(3)
 out+=['''<p>先把模型本体写清楚：方向键消息传递、均值读出、回归头三部分，分别对应“局部更新”“汇总成分子”“输出一个数”。完整实现见common.py的make_model。</p>''' ,code('''from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN
from chemprop.models import MPNN
model = MPNN(
    BondMessagePassing(d_h=128, depth=3, dropout=0.05),
    MeanAggregation(),
    RegressionFFN(n_tasks=1, input_dim=128, hidden_dim=128,
                  n_layers=2, dropout=0.05))'''),'''<p>代码实践里最值得防的坑，不是“网络不够深”，而是<strong>预测到底在哪个量纲上</strong>。pIC50在训练子集内标准化，均值和标准差只能从该子集计算。验证和测试预测再用同一组统计量还原，不能重新计算，也不能反变换两次。</p>
<p>下面是完整脚本中实际使用的核心逻辑。为了让张量流向更直观，本篇用PyTorch训练循环调用真实Chemprop模型；它不是官方CLI默认训练流程的逐项复刻。训练损失是标准化目标上的MSE，选轮次则看原尺度验证MAE。</p>''',code('''# scripts/e19_benchmark.py：同一子集的统计量只计算一次
mean = float(y[train_idx].mean())
scale = float(y[train_idx].std(ddof=0))

# 模型没有挂输出反变换层，输出仍是标准化尺度
pred_norm = model(batch.bmg, batch.V_d, batch.X_d)
target_norm = (batch.Y - mean) / scale
loss = ((pred_norm - target_norm) ** 2).mean()
loss.backward()
optimizer.step()

# 评估：eval模式、无梯度，再手动还原一次
with torch.no_grad():
    model.eval()
    pred_norm = model(batch.bmg, batch.V_d, batch.X_d)
    pred_pIC50 = pred_norm.cpu().numpy() * scale + mean'''),'''<p>上面是节选，完整脚本还包含zero_grad、学习率调度、梯度裁剪、随机种子、验证与检查点恢复。若使用Chemprop的自动输出变换，应按版本文档配置UnscaleTransform，就不要再手动乘回一次。<a href="https://chemprop.readthedocs.io/en/latest/tutorial/python/scaling.html">官方缩放说明</a>。</p>''',code('''# 完整可复跑入口：在配套目录根下执行
python scripts/e18_prepare.py
python scripts/e19_benchmark.py

# 先检查小数据通路；配置仍读取同一个config.json
python scripts/e19_benchmark.py --sizes 500 --seeds 42'''),'''<p>训练循环每轮查看验证集，不查看测试集；停止后恢复验证集最优权重，再统一导出测试预测。3个种子全部报告，不挑其中“最好看”的一轮或一个种子当作主成绩。<a href="https://chemprop.readthedocs.io/en/latest/training.html">Chemprop模型构建与训练接口参考</a>。</p>''',fig(3,'图03_验证轨迹.png','7148分子拟合池的3条验证轨迹。星号标记每次最终保留的最佳验证轮次；曲线不代表测试集调参。')]
 h(4)
 rows=[]
 for model,label in [('RF','Morgan＋RF'),('GNN','D-MPNN')]:
  r=lookup[full,model];rows.append([label]+[f"{r[k+'_mean']:.4f} ± {r[k+'_std']:.4f}" for k in ['MAE','RMSE','R2','spearman']])
 out+=[table(['模型','MAE↓','RMSE↓','R²↑','Spearman↑'],rows),f'''<p>读表时先分清任务：MAE与RMSE越低越好，R²与Spearman越高越好；Spearman检查排序单调性，不保证具体活性数值准确。MAE用pIC50的对数单位计量，也不是错误率。</p>
<p>上表是3次运行的<strong>均值±样本标准差</strong>。本轮拟合池为{full}个分子，GNN有86273个可训练参数；{winner}的平均MAE更低。差值是{abs(delta):.4f}个pIC50单位，不是“准确率差了{abs(delta)*100:.1f}%”。</p>''',fig(4,'图04_同测试集散点.png','seed42的逐分子预测，两边用同一测试集合与坐标范围。只展示一个种子便于读图，主成绩仍是上表的3次汇总。'),'''<p>为什么还要看每个种子？平均值可能掩盖不稳定。把同一训练子集与种子下的两模型相减，定义ΔMAE = GNN − RF，负值代表该次GNN误差更低：</p>''',table(['种子','ΔMAE：GNN − RF','GNN最佳轮次／已运行轮次'],[[seed,f'{value:+.4f}',str(int(full_gnn.loc[int(seed),'best_epoch']))+'／'+str(int(full_gnn.loc[int(seed),'epochs_executed']))] for seed,value in summary['full_seed_paired_delta_GNN_minus_RF'].items()]),'''<p>这里也报告最佳轮次，是为了检查训练预算是否成为限制：若最优验证轮次靠近50轮上限，就不能宣称模型已经充分收敛，更不能直接把劣势归因为“数据量不够”。学习率、训练轮次与架构仍是未分离的因素。</p>
<p>这里也要避免另一种误读：第11篇RF的骨架簇bootstrap区间，是固定该模型与测试集合的条件区间，不能用“GNN分数落在里面”判断两模型差距是否显眼。<strong>要研究两模型差值的不确定性，应直接对同一测试分子上的成对误差差值建区间。</strong>本轮补充复核已直接构造骨架簇配对差值区间，见下文；三个种子的范围仍不能代替该区间，更不能代表更换训练集或测试划分后的风险。</p>
<p><strong>标准差不等于模型已经通过了泛化检验。</strong>我们没有更换外层骨架划分，所以这些重复主要反映初始化、数据顺序和子集抽样的波动，不包含换一批新骨架后的变化；3个种子也不足以精确估计尾部风险。</p>
<p>第11篇的旧快照中，指纹＋RF的MAE为0.6015，408601参数、训练50轮的Chemprop为0.6394。这两个旧结果不能直接放进本轮表格争名次：旧模型使用8698个训练分子，新协议留了1550个作验证，而且网络配置与训练实现也变了。<strong>本轮对照必须使用本轮重新训练的RF。</strong>差异可能来自模型、特征、优化和验证方案；图特征包含手性编码，而本轮RF位指纹未启用手性，输入信息也不完全对称。不能把整段差值都叫作“图表示的纯收益”。</p>''']
 out+=['''<h3 id="attachment-protocol-audit">附件里的40万参数成绩，为什么与本轮不同？</h3>
<p>附件提供了另一组真实的300维D-MPNN快照。我们逐一读取16份NPZ预测，检查测试索引、真值和指标，均与g1对应汇总记录吻合。<strong>它值得作为补充证据，但不是本轮128维实验的替代结果，也不是本轮重新训练出来的结果。</strong></p>''',table(['协议','拟合／验证／测试','GNN与选择规则','RF／GNN的平均MAE'],[
['A：本文主实验','7148／1550／1770','128维；验证集选轮次；3种子',f"{rf['MAE_mean']:.4f}／{gnn['MAE_mean']:.4f}"],
['B：附件骨架快照','8698／不另设验证／1770','300维；固定50轮取末轮；5种子',f"{extra['rf']['mean']:.4f}／{extra['gnn']['mean']:.4f}"]
]),'''<p>两套协议共用骨架测试分子，但训练数量、隐藏维、批大小、训练循环、种子集合和模型选择规则均不完全相同。<strong>不能把协议A到B的变化，单独归因于隐藏维或数据量。</strong>“RF可学参数为0”和“唯一变量只是分子表示”也都不成立：树结构与叶节点输出同样由数据学习。</p>''',supfig('A','图09_双协议与配对差值.png','上方把各次运行逐点展示，避免只看均值；下方对同一测试分子的成对误差差值做骨架簇bootstrap。协议B是附件预测复核，不是本轮重训。'),table(['固定预测集合','平均差值：GNN−RF','95%骨架簇条件区间'],[
['A：3种子误差均值',f"{paired['main_128dim_3seeds']['mean_delta_GNN_minus_RF']:+.4f}",f"[{paired['main_128dim_3seeds']['cluster_bootstrap95'][0]:+.4f}, {paired['main_128dim_3seeds']['cluster_bootstrap95'][1]:+.4f}]"],
['B：5种子误差均值',f"{paired['attachment_300dim_5seeds']['mean_delta_GNN_minus_RF']:+.4f}",f"[{paired['attachment_300dim_5seeds']['cluster_bootstrap95'][0]:+.4f}, {paired['attachment_300dim_5seeds']['cluster_bootstrap95'][1]:+.4f}]"]
]),'''<p>具体做法是：先对每个测试分子的各次绝对误差差值取均值，再以测试集的771个骨架簇为单位，同步抽取两模型的误差，重复2000次。这里平均的是误差，不是先把预测做集成。两套协议的条件区间均在0以上，支持“在当前固定预测与测试骨架上，GNN平均误差更大”的观察；<strong>区间不包含重新训练、重新划分或新项目的不确定性，也不是模型家族普遍优劣的证明。</strong></p>
<p>附件还提供了g1的随机划分3种子预测：这些文件也已逐一核对。不过另一个g4批次使用相同种子时产生了不同记录，例如seed42的GNN MAE分别为0.5421和0.5490。我们没有混合两个批次，也不会在缺少完整运行环境记录时，仅凭数字不同就归因为GPU不确定性。</p>
<p>附件的5种子标准差比值约为4.9，描述的是这组运行的离散程度，<strong>不是“方差5倍”或实验测量噪声的因果解释</strong>。计时同样要核对：该批GNN记录中存在3213.3秒的长运行，不能用“每次约7分钟”覆盖全部情况，或事后无说明地删掉它。</p>''',code('''# scripts/e23_attachment_audit.py：检验差值，保留同一分子配对
per_seed_delta = abs(pred_gnn - y_true) - abs(pred_rf - y_true)
# 多种子时先在同一分子上取平均，再整组抽取测试骨架
# 每次同时抽取两模型对应误差，计算加权平均差值
# 完整实现见cluster_interval；输出为条件区间，不是外部验证。''')]

 h(5)
 curve_rows=[]
 for n in sizes:curve_rows.append([int(n),f"{lookup[n,'RF']['MAE_mean']:.4f}",f"{lookup[n,'GNN']['MAE_mean']:.4f}",f"{lookup[n,'GNN']['MAE_mean']-lookup[n,'RF']['MAE_mean']:+.4f}"])
 out+=['''<p>如果“GNN只是数据不足”是一个值得检验的解释，至少要看看数据增加时，它的误差怎样变化。每个种子先打乱拟合池，取前500、前2000、前5000、最后全部7148个分子；这样小集合包含在大集合里，能减少不同规模完全换了一批分子的干扰。这里是按分子随机取样，不是按骨架整组抽样；样本量增加时，骨架覆盖也会变化，所以仍不能把全部变化解释成纯粹的数量效应。</p>''',code('''# 两种模型使用同一个train_idx；验证/测试始终固定
order = np.random.default_rng(seed).permutation(fit_pool)
for n in [500, 2000, 5000, len(fit_pool)]:
    train_idx = order[:n]
    # 分别训练RF和GNN；不得把val/test拼回去凑样本量
    # 完整实现见scripts/e19_benchmark.py'''),fig(5,'图05_学习曲线.png','实跑24次拟合，每个点汇总3个种子；阴影是均值±样本标准差，不是置信区间。'),table(['实际拟合分子','RF平均MAE','GNN平均MAE','ΔMAE'],curve_rows),f'''<p><strong>{curve_statement}</strong>曲线只覆盖这个靶点的当前数据范围，不能给出“5000个分子必须用RF、十万个分子一定用GNN”的通用界线。增加数据没有保证某种架构一定反超；固定50轮与预设学习率也可能限制了优化，尚不能分离所有原因。</p>''','<h3>学习曲线的四种读法</h3>',table(['看到的曲线','可以怎样写','不能直接推出什么'],[
['GNN在各采样点的平均MAE都更高','当前数据范围与预设配置中，未观察到GNN的平均MAE优势','GNN不适合这个任务，或图表示没有价值'],
['采样点之间发生均值排名变化','报告哪些采样规模下排名不同，并同时展示种子波动','连线交点就是精确数据门槛，或超过门槛必然获胜'],
['GNN误差下降更快，但仍未追平','在已观察区间内下降更快，尚未追平','继续加数据就一定反超'],
['一个模型看似趋平，另一个仍下降','当前训练方案与采样范围下，曲线形态不同','某类模型已到理论上限，或更深就会更好']
]),'''<p>附件另有502、1001、2000、5086、8006、8698分子的曲线，整组抽取训练骨架，但每档只有seed42；中间档仅附JSON汇总，没有逐分子预测供独立复算。它可作为另一个采样协议的探索记录，不能把全量档的5种子波动直接搬成所有档位的误差带。本篇主图继续采用实际拟合池中的4档×3种子结果。无论按分子抽样还是按骨架组抽样，数量增长都可能同时改变化学覆盖，不能宣称自动消除了相似度混淆。</p>
<p><strong>本次对应第一种观察：</strong>四个采样规模上，GNN的平均MAE均高于RF。我们可以保留这个负结果，但不能把它升级成模型家族的终局判断。种子波动应与均值一起报告，区间是否重叠也不能单独替代配对比较。</p>''',fig(6,'图06_训练成本.png','同一台机器实测的训练流程墙钟时间。RF仅计拟合、使用4个并行作业，GNN计拟合及逐轮验证、使用1个Torch线程；时间不等于同硬件同线程下的算法复杂度比较。'),f'''<p>在全拟合池上，RF平均用时<strong>{rf['fit_seconds_mean']:.1f}秒</strong>，GNN平均<strong>{gnn['fit_seconds_mean']:.1f}秒</strong>，本机这组资源配置下约为{ratio:.1f}倍。RF只计拟合，GNN计拟合及选择轮次所需的逐轮验证；排除了共同预处理和最终测试推断；部分种子任务并行执行，墙钟时间还受资源争用影响。它可用于这台电脑的运行规划，不能当作CPU与GPU的通用速度排名。</p>
<p>再提醒一个常见误区：“RF参数少，所以一定更不容易过拟合”并不成立。树数量、叶节点复杂度、特征、训练集大小都会影响RF；神经网络的权重数和森林的树节点数也不是同一把复杂度尺。选模型要看验证结果和资源，而不是只看“传统”还是“深度”。</p>''']
 out+=['''<h3>模型升级，最后要算三笔账</h3>''',table(['预算','这次应该看什么','不能用什么代替'],[
['数据预算','学习曲线、化学系列覆盖与标签质量；固定指纹虽无需拟合提取规则，RF预测器仍要学习','拿神经网络权重数÷分子数，直接判断数据够不够'],
['算力预算','实际训练流程用时、硬件、线程和验证开销；收益能否覆盖项目成本','把不同资源与选择规则下的墙钟比值当作通用速度排名'],
['解释预算','能否从预测回查特征、环境与失败案例；归因方法还有哪些限制','把指纹位映射成唯一药效团，或把GNN热力图当作因果证据']
]),'''<p>这三笔账没有统一权重。若项目只需要尽快建立可靠排序，维护成本和复跑速度可能很重要；若要研究新的表示方案，则需要为更充分的验证与诊断预留预算。<strong>是否值得升级，要由具体用途和证据决定。</strong></p>''']

 h(6)
 out+=['''<p>成绩单告诉我们谁的平均误差更小，却没有告诉我们模型对哪些输入变化更敏感。现在给模型做一次“输入体检”：<strong>同一个分子，只改动一个原子的输入特征，预测会怎样变？</strong></p>
<p>这能帮助我们读懂模型行为，但不能单凭一张热图判断“它已经学懂化学”或“只是背题”。本次Chemprop的D-MPNN没有注意力层，所以这里采用<strong>节点特征遮蔽敏感性</strong>：把一个原子的特征向量整体置零，保留全部原子、化学键和连接，观察输出变化。</p>
<h3>先说明这张图对应哪一个模型</h3>''',f'''<p>本节主案例使用<strong>协议A、7148分子拟合、128隐藏维、seed42</strong>的模型。厄洛替尼完整输入的预测pIC50为{readout['base_prediction']:.4f}；它位于拟合池中，只作解释演示。后面附件300维模型的7.7220属于另一套协议，两个预测不能混用。</p>''',code('''# e25用masked-full，便于直接读“遮蔽之后升高还是降低”
# e20保存的是full-masked；e25显式反转符号，两者不混用。
masked = copy.deepcopy(batch_graph)
masked.V[atom_index] = 0
with torch.no_grad():
    full_pred = model(batch_graph).item() * scale + mean
    masked_pred = model(masked).item() * scale + mean
delta = masked_pred - full_pred'''),'''<p>对照的RF来自同一次主实验。取其<strong>全局重要性最高的20个指纹位</strong>，通过RDKit的bitInfo映射检查哪些环境出现在厄洛替尼中，给对应环境的中心原子计票。它与GNN的单分子遮蔽不是同一种指标，放在一起是为了展示两条检查线索，不是比较谁“更懂药效团”。</p>''',code('''# scripts/e25_case_readout.py：环境中心计票
extra = rdFingerprintGenerator.AdditionalOutput()
extra.AllocateBitInfoMap()
generator.GetFingerprint(mol, additionalOutput=extra)
bit_info = extra.GetBitInfoMap()
for bit in global_top20_bits:
    for center, radius in bit_info.get(int(bit), []):
        center_votes[center] += 1
# 这只统计中心，不等于该环境中所有原子的贡献。
# 完整脚本同时保存环境覆盖原子列表，供逐位回查。'''),fig(7,'图07_敏感性与指纹位.png','左为协议A逐原子特征遮蔽后的预测变化，右为绝对变化最大的10个原子，旁注RF全局Top20位的环境中心计数；编号从0开始。红色表示遮蔽后预测升高，蓝色表示降低。RF票数不是局部活性贡献。'),'<h3>读图发现一：敏感性分布在多个节点上</h3>',f'''<p>29个重原子中，有<strong>{readout['n_abs_delta_gt_threshold']}个</strong>的遮蔽绝对变化超过0.02个pIC50单位；最大变化约<strong>{readout['max_abs_delta']:.3f}</strong>，对应中心NH氮原子N{readout['most_sensitive_atom']}。这说明这一次输出对多个节点的输入编码有响应，而不是只能由一处输入改变。</p>
<p>0.02只是这里用于描述的阈值，不是统计显著性或药理阈值。把特征置零也不等于删掉原子、合成类似物，不能把“变化0.575”读成“这个氮贡献了0.575活性”。</p>''','<h3>读图发现二：两份名单部分重合，但不能直接作机制比较</h3>',f'''<p>RF映射有<strong>{len(readout['rf_center_voted_atoms'])}个</strong>中心原子获得票数；GNN最敏感的10个原子中，有<strong>{len(readout['top10_and_center_votes_overlap'])}个</strong>也在这份中心名单中。这个结果既不是“完全读同一块”，也不是“几乎不重叠”。</p>
<p>更关键的是，重合度同时受全局Top20位筛选、哈希碰撞、中心/环境映射口径以及遮蔽方法影响。全局不纯度重要性不是这个分子的局部归因；中心票也不是整个子结构的全部原子。名单可用于提出检查问题，不能证明模型学到了不同药效机制。</p>''','<h3>读图发现三：它给出线索，没有解释模型输赢的原因</h3>', '''<p>可以继续问：更换读出方式会不会改变敏感性？有意义的结构替换能否被模型正确响应？但要回答这些问题，需要单独设计消融或使用有实验标签的匹配分子对。<strong>本图不能把主成绩差距归因为“GNN容量过大”“指纹提供强正则”或“模型没学会化学”。</strong>未见分子上的表现，仍要回到测试集、学习曲线和独立验证。</p>
<p>零向量可能超出训练分布，所以这里的颜色是输入扰动诊断，不是实验结论。先把这条边界守住，热力图才不会比模型本身更会讲故事。</p>''','<h3 id="scaling-audit">实验加餐：同一权重，为什么会从9.82变成7.72？</h3>',f'''<p>附件归因图报告厄洛替尼的预测pIC50为9.8224。复核代码发现：模型训练时使用全部8698个训练分子拟合目标缩放器，但归因时只取前200个训练分子重新计算均值与标准差。两者并非同一把尺。</p>
<p>本次用同一300维权重和正确训练统计量重新推断全部1770个测试分子，与附件保存预测的最大差异为{attachment['model_reload_max_abs_difference']:.2g}，确认了还原路径。随后只修正归因时的缩放统计量，厄洛替尼的预测由<strong>{erlo['original_reported_prediction']:.4f}修正为{erlo['correct_scale_prediction']:.4f}</strong>。</p>''',supfig('B','图10_缩放还原审计.png','同一附件300维模型的预测，还原尺度不同就会改变数值。修正使用训练时的统计量；这不是重训模型，也不是一次模型性能升级。'),table(['缩放统计量','错误归因路径：前200个','训练时实际使用：8698个'],[
['均值',f"{correction['wrong_mean']:.4f}",f"{correction['correct_mean']:.4f}"],
['标准差',f"{correction['wrong_scale']:.4f}",f"{correction['correct_scale']:.4f}"]
]),'''<p>遮蔽差值也会被标准差缩放。因此，附件图中的0.21等原尺度敏感性数值不能直接沿用。我们已将修正后的三分子结果保存到复核JSON，同时将本篇协议A的归因图与结果分别报告，两者不混用。</p>
<p>另一个问题来自分子身份：附件示例SMILES没有在原始字符串列中匹配成功，但规范化后与数据第3880行相同，位于训练集中。<strong>“字符串不同”不等于“新分子”。</strong>本篇保留它作解释演示，不把它称作训练集外验证。</p>''',code('''# 正确做法：训练时保存缩放器或其参数，解释时直接复用
checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
with torch.no_grad():
    pred_norm = model(batch_graph)
pred_pIC50 = pred_norm * checkpoint["scale"] + checkpoint["mean"]
# 不要用“前200个”或待预测分子重新估计mean/scale。
# 附件旧权重的训练统计量由e23按原训练索引重建，并用保存预测核对。'''),'<h3>再回到未参与拟合的分子：局部解释与预测效果分开看</h3>',fig(8,'图08_成对误差.png','seed42固定测试集上的成对误差；右图负值表示该分子上GNN误差较小。训练池内分子的热力图不能代替这张测试成绩图。'),f'''<p>在这批测试分子上，GNN有{summary['seed42_fraction_GNN_lower_abs_error']:.1%}的分子绝对误差小于RF，但总体平均MAE仍更大。可见“某些分子预测得更好”“总体更准”和“解释看起来合理”是三个不同判断。下面的案例按固定规则选取，只用于事后诊断，不作为独立验证：</p>''',table(['案例选择规则','实验pIC50','RF预测','GNN预测'],[[r['rule'],f"{r['true']:.3f}",f"{r['RF']:.3f}",f"{r['GNN']:.3f}"] for r in summary['cases']])]
 h(7)
 out+=['''<h3>负结果也是结果，但要把主语写完整</h3>
<p>本次没有观察到GNN的平均MAE优势。值得保留的结论是：<strong>在这份EGFR快照、固定骨架测试集、128维D-MPNN与当前训练方案下，指纹＋RF是表现更好的基线。</strong>不是“GNN不行”，也不是“数据一定不够”。</p>
<p>我们报告全部三个种子、四个规模和训练轨迹，没有只留下最好看的一次。如果下一步要调整深度、读出、手性配置或训练轮数，应另立实验版本，记录变更并在验证集上作选择；不能在同一测试集上反复试到翻盘，再只发布翻盘的结果。已有历史测试集已被多次分析，更广泛的结论还需要独立数据验证。</p>
<p>本篇保存了正式配置、数据哈希和训练记录；这与公开、带时间戳的事前注册不是同一件事，也不构成“从未见过测试结果”的保证。</p>
<h3>回到项目：看到什么证据，采取什么动作？</h3>''']
 out+=['''<p>把精度、数据与运行/解释成本放在一起，才能回答“我现在该不该上GNN”。下面按实际场景给出建议，每条都写明哪些是本篇证据、哪些仍需另做验证：</p>''',table(['你的场景','建议采取的动作','依据与适用边界'],[
['单靶点，几千到万级数据','先建立指纹＋树模型基线，再在同一评估框架下试GNN','本篇两套预设都未显示GNN平均MAE优势；不是所有单靶点的通用排名'],
['训练数据还在增加','重复学习曲线，同时记录化学覆盖、波动与成本','本篇范围内没有均值排名反转，不据此猜测未来交叉点'],
['多靶点、预训练或迁移任务','把图模型与其他表示列为候选，另设任务匹配的对照','本篇没有测试迁移收益，不能称其为已经证明的GNN主场'],
['需要三维姿势或几何信息','先确认输入是否包含相应信息，再为结构方法单独验证','本篇两个2D模型均未输入坐标；换成3D也不自动证明活性排序可靠'],
['需要向化学家解释预测','同时提供指纹位环境回查、输入敏感性与失败案例','可追溯性有价值，但中心票和热力图都不等于药效机制'],
['需要区分手性或立体选择性','检查SMILES和特征的手性编码，再用相关实验数据验证','能编码两种构型，不等于已经学对两者的活性差异'],
['算力、时间或维护预算有限','按验证表现、稳定性与运行成本，保留足够好的简单基线','RF和GNN的容量形式与计时口径不同，不能用“参数0对40万”定夺']
]),'''<p>在流水线中，RF与GNN是<strong>活性预测环节的两种可选实现</strong>。结构准备与对接提供可审阅的姿势假设，活性模型提供历史数据支持的优先级，ADMET相关模型评估其他属性；它们回答不同问题。把几项输出放在一起，不代表可以直接混合单位，更不能把实验真值写回候选排序。</p>
<blockquote><p><strong>模型选择不是信仰问题，是证据与预算问题。</strong>先把比较规则说清，再看数据范围内的增益、波动与代价；值得升级时升级，没有证据时继续保留基线。</p></blockquote>
<p>本篇还有明确的未完成范围：只研究单靶点和一种外层骨架划分，没有超参数搜索、独立项目验证、多靶点迁移或3D图模型。测试集合在第11篇已被分析，本篇属于回顾性基准实验；若后续根据这张成绩单继续挑配置，就应另找新的独立验证数据。</p>
<h3>向后回接第8篇：升级的是MVP中的一个模块</h3>
<p>回到第8篇跑通的流程，本篇的交付不是“所有模型都换成GNN”，而是一种可复用的升级方法：在固定输入和评价规则下，比较候选预测器，保留模型版本、训练来源与验证记录，再决定是否替换活性预测模块。本次配置没有显示GNN的平均MAE优势，继续保留RF基线，就是有实验依据的选择。</p>
<p>蛋白结构与对接提供结合姿势假设，活性预测安排候选优先级，ADMET与合成可行性评估其他约束；这些模块可以并行、交叉检查并在新实验后更新。它们的输出不能互相代替证据，也不能直接把不同单位与未经校准的分数相加。</p>
<h3>向前连接第13篇：从评分分子，到提出新分子</h3>
<p>按当前学习规划，下一篇进入分子生成：从“给一个已有分子打分”转向“提出哪些新结构值得检查”。图表示在这里提供理解原子、键与结构约束的基础，但生成模型不一定都采用分子图，也可能使用SMILES等表示。</p>
<p><strong>生成器给出的高分，首先只是评分模型的高分。</strong>如果候选离训练分布越来越远，生成器还可能利用评分模型的盲区。第11篇的划分与适用范围、本篇的对照规则和误差诊断，都会成为下一篇筛查生成候选的工具；生成结果仍须经过化学有效性、属性、可合成性与实验等层次的检查。</p>
<p>ADMET、逆合成与干湿闭环仍是后续规划主题，具体展开顺序再按系列进度确定。整条学习路径最终回到同一个目标：<strong>让新证据持续回流，逐步把第8篇“能跑通”的MVP，升级成每个关键环节都有依据的流程。</strong></p>''']
 out+=['<h2 id="reproduce">配套目录、复跑命令与资料来源</h2>','''<pre class="directory">AI制药第12篇_图神经网络/
├── AI制药第12篇_公众号文章.html   # 通过相对路径引用配图
├── 优化提纲.md                    # 参考规划的取舍与优化
├── config.json                    # 正式模型配置与种子
├── data/                          # 固定分子、指纹、划分与历史摘要
├── scripts/
│   ├── common.py                  # 路径、固定协议、模型构建
│   ├── e18_prepare.py             # 数据/特征/版本审计
│   ├── e19_benchmark.py           # 24次拟合，可断点续跑
│   ├── e20_sensitivity.py         # 节点遮蔽与RF位映射
│   ├── e21_figures.py             # 汇总结果、生成8张图
│   ├── e22_article.py             # 结果驱动生成正文与内嵌预览
│   ├── e23_attachment_audit.py    # 附件预测复核、配对区间与缩放纠错
│   ├── e24_supplement_figures.py  # 重绘补充图9、10
│   ├── e25_case_readout.py        # 主案例敏感性与RF中心映射读图
│   ├── make_course_map.py         # 重绘课程定位导读图
│   └── verify_article.py          # 交付一致性检查
├── attachment_evidence/           # 附件原稿、脚本、预测及一份权重快照
├── results/                       # 主实验指标与附件复核结果
├── models/gnn_full_seed42.pt       # 主模型权重与缩放统计量
├── pilot_300dim/                  # 通路检查，不纳入主结果
├── figures/                       # 正文图1—8与独立概念封面
├── requirements.txt               # 本次依赖版本
└── README.md                      # 环境、运行方法与局限
</pre>''',code('''# 使用已安装requirements.txt依赖的Python；在配套目录根运行
python scripts/e18_prepare.py
python scripts/e19_benchmark.py
python scripts/e20_sensitivity.py
python scripts/e21_figures.py
python scripts/e23_attachment_audit.py
python scripts/e24_supplement_figures.py
python scripts/e25_case_readout.py
python scripts/make_course_map.py
python scripts/e22_article.py
python scripts/verify_article.py'''),'''<p>e19会跳过已完成且配置哈希一致的运行；改配置后应使用新结果目录，避免把不同实验混在一起。所有输入快照、训练子集索引、数据哈希、预测与轨迹均随包提供。</p>
<p>资料来源：<a href="https://pubmed.ncbi.nlm.nih.gov/31361484/">Yang等：Analyzing Learned Molecular Representations for Property Prediction（2019）</a>；<a href="https://chemprop.readthedocs.io/en/latest/training.html">Chemprop训练文档</a>；<a href="https://chemprop.readthedocs.io/en/latest/tutorial/python/scaling.html">Chemprop输入输出缩放文档</a>。接口以本机2.3.1为准，文档可能更新；固定版本和输入快照见配套包。</p>
<p class="meta">本篇为公开分子数据上的教学实验。协议A来自本轮24次拟合；协议B来自附件既有预测的复算与权重重载，未在本轮重新训练；配对重抽样和缩放纠错为本轮新增诊断。0.6015与0.6394为前期工程快照，均已标明来源。ChEMBL数据库发布号未随原始快照记录，本包支持固定输入复跑，不承诺从在线最新版重新抓取获得相同样本。预测不替代实验确认。</p>
<section id="series-index"><h2>系列文章索引｜零基础入门AI制药</h2><p>先建立全局地图，再沿MVP的各个模块深入。下面按本系列课程导航排列；标题为便于查找的主题简称。</p>
<h3>第一阶段：认识问题、准备工具、跑通流程</h3>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第1篇</span> · <a href="https://mp.weixin.qq.com/s/z-NGKlnCenUDgTdvwxs0Eg">起始篇：为什么需要AI制药</a><br><span style="font-size:13px;color:#647571;">从新冠mRNA到计算设计，建立学习地图</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第2篇</span> · <a href="https://mp.weixin.qq.com/s/sL5O-O6CanABKy9-De1HhA">第1课：和AI读懂一个分子</a><br><span style="font-size:13px;color:#647571;">认识SMILES、分子指纹与基础分子语言</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第3篇</span> · <a href="https://mp.weixin.qq.com/s/-HGDZi7s8nCA-c-jU0KcPg">第2课：讲透四层生化语义</a><br><span style="font-size:13px;color:#647571;">理解RDKit输出背后的化学与生物含义</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第4篇</span> · <a href="https://mp.weixin.qq.com/s/lG1t9PLBeOuWby7I1PkuDA">政策篇：创新药的政策命题</a><br><span style="font-size:13px;color:#647571;">理解研发目标、AI的作用与边界</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第5篇</span> · <a href="https://mp.weixin.qq.com/s/6k62fZ5aAaNOzxPr84aD0w">工具篇：从RDKit到DeepChem</a><br><span style="font-size:13px;color:#647571;">搭建AI制药工具与实践环境</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第6篇</span> · 实践篇·上：建立全流程地图<br><span style="font-size:13px;color:#647571;">从疾病、靶点到候选分子，认识各环节的接口</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第7篇</span> · 行业篇：理解产业全景<br><span style="font-size:13px;color:#647571;">认识市场、技术路线与商业模式</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第8篇</span> · 实践篇·实操：跑通端到端MVP<br><span style="font-size:13px;color:#647571;">亲手执行流程，检查数据与工程接口</span></p>
<h3>第二阶段：拆开MVP，逐个检查技术模块</h3>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第9篇</span> · 蛋白质篇：结构、可信度与结合口袋<br><span style="font-size:13px;color:#647571;">获取结构，区分结构可信度与位点判断</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第10篇</span> · <a href="https://mp.weixin.qq.com/s/gb-3thG7LOs0u94MrXFrJw">三维空间与分子对接：从Vina到DiffDock</a><br><span style="font-size:13px;color:#647571;">区分姿势预测与活性排序，检查评分边界</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第11篇</span> · 分子表征与活性预测：指纹＋随机森林<br><span style="font-size:13px;color:#647571;">建立QSAR基线，用随机/骨架划分审视泛化</span></p>
<p style="margin:12px 0;"><span style="color:#0F755B;font-weight:700;">第12篇</span> · <strong>图神经网络与分子表征升级（本篇）</strong><br><span style="font-size:13px;color:#647571;">比较指纹与图模型的精度、波动、解释和成本</span></p>
<p class="meta">带链接的篇目可点击阅读；其余篇目可在公众号历史消息中按主题查找。</p><h3>接下来：从评价分子，到生成与验证分子</h3><p>第13篇计划讨论分子生成与约束优化；ADMET、逆合成与干湿闭环为后续规划主题，具体发布顺序以系列更新为准。</p><blockquote><p>推荐阅读路线：第8篇看全流程 → 第9—10篇理解结构与对接 → 第11—12篇检查活性预测 → 后续把评价规则带入分子生成。</p></blockquote></section>
</main></body></html>''']
 # 正文的外部参考文献保留文字，仅公众号链接可点击；本地预览保留目录锚点。
 def clean_link(match, keep_anchors=True):
  url=html.unescape(match.group(1))
  if (keep_anchors and url.startswith('#')) or (urlsplit(url).scheme=='https' and urlsplit(url).hostname=='mp.weixin.qq.com'):
   return match.group(0)
  return '<span>'+match.group(2)+'</span>'
 text='\n'.join(out)
 text=re.sub(r'<a\b[^>]*href="([^"]+)"[^>]*>(.*?)</a>',clean_link,text,flags=re.S)
 article=ROOT/'AI制药第12篇_公众号文章.html';article.write_text(text)
 publish=re.sub(r'<a\b[^>]*href="([^"]+)"[^>]*>(.*?)</a>',lambda m:clean_link(m,False),text,flags=re.S)

 for ref in re.findall(r'<img[^>]+src="([^"]+)"',text):text=text.replace('src="'+ref+'"','src="data:image/png;base64,'+base64.b64encode((ROOT/ref).read_bytes()).decode()+'"')
 (ROOT.parent/'AI制药第12篇_完整预览.html').write_text(text)
 for ref in re.findall(r'<img[^>]+src="([^"]+)"',publish):publish=publish.replace('src="'+ref+'"','src="data:image/png;base64,'+base64.b64encode((ROOT/ref).read_bytes()).decode()+'"')
 (ROOT.parent/'AI制药第12篇_公众号发布版.html').write_text(publish)
 # 摘要不预设赢家；紧随其后的结论从实跑结果生成。
 abstract='把2048位指纹换成分子图，预测更准吗？本文沿用EGFR骨架划分，完成24次模型拟合，对比随机森林与Chemprop的精度、波动和成本。结合学习曲线与原子特征遮蔽，检验模型升级的收益与局限，解释为何更复杂未必更好。附脚本、结果和十一张配图。'
 (ROOT/'发布资料.md').write_text(f'# 推荐标题\n\nAI制药第12篇：{headline}\n\n## 120字摘要\n\n{abstract}\n\n## 本轮结论\n\n全拟合池3种子平均MAE：RF {rf["MAE_mean"]:.4f}，GNN {gnn["MAE_mean"]:.4f}。{curve_statement}\n\n## 发布检查\n\n正文代码为黑底。内嵌预览可单独打开；公众号后台需上传11张独立图片，并复查手机表格、目录链接与代码换行。未自动发布或上传仓库。\n\n## 封面\n\nfigures/公众号封面_思路与摘要版.png；内置imagegen生成，提示词见封面生成提示词_思路摘要版.txt。封面为概念示意，正文8张主图、2张补充实验图与1张课程定位图由脚本绘制。\n')
 print('wrote article, embedded preview, publication notes')
if __name__=='__main__':main()
