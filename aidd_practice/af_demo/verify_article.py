#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：交付校验：文章 HTML
# 使用注解：检查文章结构与引用等，不替代数值科学审查；执行前确认目标路径。
# 完整对应关系与运行顺序：../README.md
"""公众号文章 HTML 交付前校验"""
import os
import re
import sys
from html.parser import HTMLParser

WORK = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
TARGET = os.path.join(WORK, "AI制药细节篇1_靶点与口袋_公众号文章.html")

VOID = {"br", "img", "hr", "meta", "link", "input", "col", "source"}
SELF_OK = {"p", "li", "td", "th", "tr", "section", "div", "h2", "h3", "pre", "table",
         "tbody", "thead", "colgroup", "strong", "code", "ol", "ul", "span"}


class Checker(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.errors = []
        self.imgs = []
        self.captions = []
        self.tables = []
        self._in_cap = False
        self._buf = ""
        self._tbl = None

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if tag == "img":
            self.imgs.append(d.get("src", ""))
        if tag == "table":
            self._tbl = {"cols": 0, "widths": []}
        if tag == "colgroup":
            pass
        if tag == "col" and self._tbl is not None:
            self._tbl["cols"] += 1
            w = d.get("style", "")
            m = re.search(r"width:\s*([\d.]+)%", w)
            self._tbl["widths"].append(float(m.group(1)) if m else None)
        if tag not in VOID:
            self.stack.append((tag, self.getpos()))

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if tag in ("table",) and self._tbl is not None:
            self.tables.append(self._tbl)
            self._tbl = None
        if not self.stack:
            self.errors.append(f"第{self.getpos()[0]}行: 多余的 </{tag}>")
            return
        if self.stack[-1][0] != tag:
            self.errors.append(
                f"第{self.getpos()[0]}行: </{tag}> 与未闭合的 <{self.stack[-1][0]}>"
                f"(第{self.stack[-1][1][0]}行) 不匹配")
            for i in range(len(self.stack) - 1, -1, -1):
                if self.stack[i][0] == tag:
                    del self.stack[i:]
                    return
            return
        self.stack.pop()

    def handle_data(self, data):
        self._buf += data
        if re.search(r"图\s*\d+\s*[｜|]", data):
            self.captions.append(data.strip())


def main():
    html = open(TARGET, encoding="utf-8").read()
    c = Checker()
    c.feed(html)

    print("=" * 62)
    print("文件:", os.path.basename(TARGET))
    print("大小: %.1f KB   行数: %d" % (len(html.encode("utf-8")) / 1024,
                                    html.count("\n") + 1))
    print("=" * 62)

    ok = True

    # 1 标签闭合
    if c.errors:
        ok = False
        print("\n[✗] 标签问题:")
        for e in c.errors[:15]:
            print("   ", e)
    elif c.stack:
        ok = False
        print("\n[✗] 未闭合标签:", [s[0] for s in c.stack])
    else:
        print("\n[✓] 标签闭合全部通过")

    # 2 图片存在
    print("\n[图片检查]  共 %d 张" % len(c.imgs))
    for s in c.imgs:
        p = os.path.join(WORK, s)
        if os.path.exists(p):
            print("   [✓] %s  (%.0f KB)" % (s, os.path.getsize(p) / 1024))
        else:
            ok = False
            print("   [✗] 缺失: %s" % s)

    # 3 图号连续
    nums = []
    for cap in c.captions:
        m = re.search(r"图\s*(\d+)\s*[｜|]", cap)
        if m:
            nums.append(int(m.group(1)))
    print("\n[图号检查]  检测到:", nums)
    if nums == list(range(1, len(nums) + 1)):
        print("   [✓] 图号连续 1–%d" % len(nums))
    else:
        ok = False
        print("   [✗] 图号不连续或有重复")

    # 4 表格宽度
    print("\n[表格检查]  共 %d 个" % len(c.tables))
    for i, t in enumerate(c.tables, 1):
        ws = [w for w in t["widths"] if w is not None]
        tot = sum(ws)
        flag = "✓" if abs(tot - 100.0) < 0.5 else "✗"
        if flag == "✗":
            ok = False
        print("   表%d: 列数=%d  宽度合计=%.1f%%  [%s]" % (i, t["cols"], tot, flag))

    # 5 微信不兼容项
    bad = []
    if re.search(r"<figure\b", html, re.I):
        bad.append("使用了 <figure>（微信不支持，请改用 <section><img><p>）")
    if re.search(r"margin:\s*-", html):
        bad.append("使用了负 margin（微信支持不稳定）")
    if re.search(r"position:\s*(fixed|sticky)", html):
        bad.append("使用了 fixed/sticky 定位")
    print("\n[微信兼容性]")
    if bad:
        ok = False
        for b in bad:
            print("   [✗]", b)
    else:
        print("   [✓] 未发现不兼容写法")

    # 6 占位符
    ph = re.findall(r"(TODO|TBD|待填|XXX|\{\{.*?\}\\})", html)
    print("\n[占位符检查]")
    if ph:
        ok = False
        print("   [✗] 发现占位符:", set(ph))
    else:
        print("   [✓] 无占位符")

    print("\n" + "=" * 62)
    print("结论:", "全部通过，可交付" if ok else "存在问题，需修复")
    print("=" * 62)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
