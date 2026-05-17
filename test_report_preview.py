#!/usr/bin/env python3
"""
生成 report.html 预览（不依赖 ViennaRNA）
"""
import os

# Mock data matching the original report structure
base_results = [
    {"rank": 1, "sequence": "AUGCGAAUUCGCGGAAACCCU", "consensus_score": 6.5, "passed": False,
     "off_target": {"risk_level": "high", "thermodynamics": {"min_dg": -2.7, "vienna_available": True,
     "seed_hits": [{"target_id": "GeneB_CYP6A2_offtarget", "dg": -2.65, "position": 3, "structure": ".(((((&)))))"}]},
     "summary": {"level_1_16bp_hits": 2, "seed_7nt_hits": 2}}},
    {"rank": 2, "sequence": "UGCGAAUUCGCGGAAACCCUU", "consensus_score": 6.3, "passed": False,
     "off_target": {"risk_level": "high", "thermodynamics": {"min_dg": -2.7, "vienna_available": True,
     "seed_hits": [{"target_id": "GeneB_CYP6A2_offtarget", "dg": -2.68, "position": 4, "structure": ".(((((&)))))"}]},
     "summary": {"level_1_16bp_hits": 2, "seed_7nt_hits": 2}}},
    {"rank": 3, "sequence": "GCGAAUUCGCGGAAACCCUUU", "consensus_score": 6.1, "passed": False,
     "off_target": {"risk_level": "high", "thermodynamics": {"min_dg": -4.2, "vienna_available": True,
     "seed_hits": [{"target_id": "GeneB_CYP6A2_offtarget", "dg": -4.18, "position": 5, "structure": ".(((((&)))))"}]},
     "summary": {"level_1_16bp_hits": 2, "seed_7nt_hits": 2}}},
    {"rank": 4, "sequence": "CGAAUUCGCGGAAACCCUUUG", "consensus_score": 5.9, "passed": False,
     "off_target": {"risk_level": "high", "thermodynamics": {"min_dg": -7.0, "vienna_available": True,
     "seed_hits": [{"target_id": "GeneB_CYP6A2_offtarget", "dg": -7.02, "position": 6, "structure": ".(((((&)))))"}]},
     "summary": {"level_1_16bp_hits": 2, "seed_7nt_hits": 2}}},
    {"rank": 5, "sequence": "GAAUUCGCGGAAACCCUUUGG", "consensus_score": 5.7, "passed": False,
     "off_target": {"risk_level": "high", "thermodynamics": {"min_dg": 0.1, "vienna_available": True},
     "summary": {"level_1_16bp_hits": 2, "seed_7nt_hits": 2}}},
]

for i in range(6, 21):
    base_results.append({
        "rank": i,
        "sequence": "AAUUCGCGGAAACCCUUUGGG"[:21],
        "consensus_score": round(5.5 - i * 0.1, 1),
        "passed": False,
        "off_target": {
            "risk_level": "high",
            "thermodynamics": {"min_dg": round(-3.0 + i * 0.1, 1), "vienna_available": True},
            "summary": {"level_1_16bp_hits": 2, "seed_7nt_hits": 2}
        }
    })

result = {
    "mode": "siRNA",
    "summary": {"total_candidates": 158, "passed_candidates": 0},
    "results": base_results
}


def _colorize_sequence(seq: str) -> str:
    color_map = {
        "A": "base-a", "a": "base-a",
        "U": "base-u", "u": "base-u",
        "T": "base-t", "t": "base-t",
        "G": "base-g", "g": "base-g",
        "C": "base-c", "c": "base-c",
    }
    return "".join(f'<span class="{color_map.get(ch, "base-n")}">{ch}</span>' for ch in seq)


def generate_html(result, output_path="report.html"):
    results = result["results"][:20]
    rows = ""
    for r in results:
        seq = r["sequence"]
        gc = round((seq.count("G") + seq.count("C")) / len(seq) * 100, 1) if seq else 0
        ot = r.get("off_target", {})
        risk = ot.get("risk_level", "N/A")
        risk_badge = {
            "low": '<span class="badge badge-success">Low</span>',
            "medium": '<span class="badge badge-warning">Medium</span>',
            "high": '<span class="badge badge-danger">High</span>',
            "N/A": '<span class="badge">N/A</span>',
        }.get(risk, f'<span class="badge">{risk}</span>')
        thermo = ot.get("thermodynamics", {})
        dg = thermo.get("min_dg", "N/A")
        dg_str = f"{dg:.1f}" if isinstance(dg, float) else str(dg)
        vienna = "✓" if thermo.get("vienna_available") else "✗"
        hits16 = ot.get("summary", {}).get("level_1_16bp_hits", 0)
        hits7 = ot.get("summary", {}).get("seed_7nt_hits", 0)
        passed = "✓ Pass" if r["passed"] else "✗ Fail"
        passed_class = "pass-yes" if r["passed"] else "pass-no"
        colored_seq = _colorize_sequence(seq)
        rows += f"""
        <tr>
            <td><span class="rank">{r['rank']}</span></td>
            <td class="mono">{colored_seq}</td>
            <td><span class="score">{r['consensus_score']:.1f}</span></td>
            <td>{gc}%</td>
            <td class="{passed_class}">{passed}</td>
            <td>{risk_badge}</td>
            <td>{hits16}</td>
            <td>{hits7}</td>
            <td><span class="dg-value">{dg_str}</span> <span class="vienna-status">{vienna}</span></td>
        </tr>
        """

    thermo_details = ""
    for r in results[:5]:
        thermo = r.get("off_target", {}).get("thermodynamics", {})
        if thermo.get("seed_hits"):
            seed = r["sequence"][1:8] if len(r["sequence"]) >= 8 else r["sequence"]
            hit_items = ""
            for hit in thermo["seed_hits"][:3]:
                hit_items += f"""
                <li>
                    <span class="hit-target">{hit['target_id']}</span>
                    <span class="hit-meta">ΔG = {hit['dg']:.2f} kcal/mol @ pos {hit['position']}</span>
                    <br><small class="hit-struct">{hit.get('structure', 'N/A')}</small>
                </li>
                """
            thermo_details += f"""
            <div class="thermo-card">
                <div class="thermo-header">
                    <span class="thermo-rank">Rank {r['rank']}</span>
                    <span class="thermo-seed">种子区 <code>{seed}</code></span>
                    <span class="thermo-dg">最小 ΔG: <strong>{thermo.get('min_dg', 'N/A')} kcal/mol</strong></span>
                </div>
                <div class="thermo-hits">
                    <p class="thermo-subtitle">热力学脱靶命中:</p>
                    <ul>{hit_items}</ul>
                </div>
            </div>
            """

    if not thermo_details:
        thermo_details = """
        <div class="empty-state">
            <div class="empty-icon">🛡️</div>
            <h3>无热力学脱靶风险</h3>
            <p>所有候选的种子区结合能均高于阈值 (-7.0 kcal/mol)，未检测到显著的热力学脱靶风险。</p>
        </div>
        """

    summary = result.get("summary", {})
    total = summary.get("total_candidates", 0)
    passed = summary.get("passed_candidates", 0)
    mode = result.get("mode", "siRNA")

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>dsRNA-Forge 设计报告</title>
    <style>
        :root {{
            --primary: #0d7377;
            --primary-dark: #0a5c5f;
            --primary-light: #14a085;
            --accent: #32e0c4;
            --bg: #f8fafc;
            --bg-card: #ffffff;
            --text: #1e293b;
            --text-secondary: #64748b;
            --border: #e2e8f0;
            --success: #10b981;
            --warning: #f59e0b;
            --danger: #ef4444;
            --radius: 12px;
            --shadow: 0 1px 3px rgba(0,0,0,0.08), 0 4px 12px rgba(0,0,0,0.05);
            --shadow-hover: 0 4px 6px rgba(0,0,0,0.1), 0 10px 24px rgba(0,0,0,0.08);
        }}

        * {{ margin: 0; padding: 0; box-sizing: border-box; }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
        }}

        /* ===== Header ===== */
        .report-header {{
            background: linear-gradient(135deg, #0f2027 0%, #203a43 50%, #2c5364 100%);
            color: white;
            padding: 48px 32px;
            position: relative;
            overflow: hidden;
        }}
        .report-header::before {{
            content: '';
            position: absolute;
            top: -50%;
            right: -10%;
            width: 500px;
            height: 500px;
            background: radial-gradient(circle, rgba(50,224,196,0.15) 0%, transparent 70%);
            border-radius: 50%;
        }}
        .header-content {{
            max-width: 1200px;
            margin: 0 auto;
            position: relative;
            z-index: 1;
        }}
        .header-badge {{
            display: inline-block;
            background: rgba(255,255,255,0.12);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255,255,255,0.15);
            padding: 6px 16px;
            border-radius: 100px;
            font-size: 12px;
            font-weight: 600;
            letter-spacing: 0.5px;
            text-transform: uppercase;
            margin-bottom: 16px;
        }}
        .header-title {{
            font-size: 36px;
            font-weight: 800;
            margin-bottom: 8px;
            letter-spacing: -0.5px;
        }}
        .header-subtitle {{
            font-size: 16px;
            opacity: 0.8;
            font-weight: 400;
        }}
        .header-meta {{
            display: flex;
            gap: 24px;
            margin-top: 24px;
            flex-wrap: wrap;
        }}
        .header-meta-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 14px;
            opacity: 0.85;
        }}

        /* ===== Stats Bar ===== */
        .stats-bar {{
            max-width: 1200px;
            margin: -32px auto 0;
            padding: 0 24px;
            position: relative;
            z-index: 2;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
        }}
        .stat-card {{
            background: var(--bg-card);
            border-radius: var(--radius);
            padding: 24px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border);
            transition: all 0.2s ease;
        }}
        .stat-card:hover {{
            transform: translateY(-2px);
            box-shadow: var(--shadow-hover);
        }}
        .stat-icon {{
            width: 40px;
            height: 40px;
            border-radius: 10px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 20px;
            margin-bottom: 12px;
        }}
        .stat-icon.mode {{ background: #e0f2fe; }}
        .stat-icon.total {{ background: #fef3c7; }}
        .stat-icon.passed {{ background: #d1fae5; }}
        .stat-icon.vienna {{ background: #ede9fe; }}
        .stat-value {{
            font-size: 28px;
            font-weight: 800;
            color: var(--text);
            margin-bottom: 4px;
        }}
        .stat-label {{
            font-size: 13px;
            color: var(--text-secondary);
            font-weight: 500;
        }}
        .stat-change {{
            font-size: 12px;
            margin-top: 8px;
            display: inline-flex;
            align-items: center;
            gap: 4px;
            padding: 2px 8px;
            border-radius: 6px;
            font-weight: 600;
        }}
        .stat-change.success {{ background: #ecfdf5; color: var(--success); }}
        .stat-change.danger {{ background: #fef2f2; color: var(--danger); }}

        /* ===== Main Content ===== */
        .main-content {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 40px 24px;
        }}

        .section {{
            background: var(--bg-card);
            border-radius: var(--radius);
            padding: 32px;
            margin-bottom: 24px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border);
        }}
        .section-header {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            margin-bottom: 24px;
            flex-wrap: wrap;
            gap: 12px;
        }}
        .section-title {{
            font-size: 20px;
            font-weight: 700;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .section-title-icon {{
            font-size: 24px;
        }}
        .section-subtitle {{
            font-size: 14px;
            color: var(--text-secondary);
            margin-top: 4px;
        }}
        .section-badge {{
            background: var(--primary);
            color: white;
            font-size: 12px;
            font-weight: 600;
            padding: 4px 12px;
            border-radius: 100px;
        }}

        /* ===== Table ===== */
        .table-wrapper {{
            overflow-x: auto;
            border-radius: var(--radius);
            border: 1px solid var(--border);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
        }}
        thead {{
            background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%);
        }}
        th {{
            padding: 14px 16px;
            text-align: left;
            font-weight: 600;
            color: var(--text);
            font-size: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            white-space: nowrap;
        }}
        td {{
            padding: 14px 16px;
            border-bottom: 1px solid var(--border);
            vertical-align: middle;
        }}
        tbody tr {{
            transition: background 0.15s ease;
        }}
        tbody tr:nth-child(even) {{
            background: #fafbfc;
        }}
        tbody tr:hover {{
            background: #f0f9ff;
        }}
        tbody tr:last-child td {{
            border-bottom: none;
        }}

        .rank {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 28px;
            height: 28px;
            border-radius: 50%;
            background: var(--primary);
            color: white;
            font-weight: 700;
            font-size: 12px;
        }}

        .mono {{
            font-family: 'SF Mono', 'Fira Code', 'JetBrains Mono', Consolas, monospace;
            font-size: 13px;
            letter-spacing: 0.5px;
            white-space: nowrap;
        }}

        /* 碱基颜色 */
        .base-a {{ color: #16a34a; font-weight: 700; }}
        .base-u, .base-t {{ color: #dc2626; font-weight: 700; }}
        .base-g {{ color: #ca8a04; font-weight: 700; }}
        .base-c {{ color: #2563eb; font-weight: 700; }}
        .base-n {{ color: var(--text-secondary); }}

        .score {{
            font-weight: 700;
            color: var(--primary);
        }}

        .dg-value {{
            font-family: 'SF Mono', monospace;
            font-weight: 600;
            color: var(--text);
        }}
        .vienna-status {{
            font-size: 12px;
            color: var(--success);
        }}

        .pass-yes {{
            color: var(--success);
            font-weight: 600;
        }}
        .pass-no {{
            color: var(--danger);
            font-weight: 500;
        }}

        .badge {{
            display: inline-flex;
            align-items: center;
            padding: 4px 10px;
            border-radius: 6px;
            font-size: 12px;
            font-weight: 600;
            background: #f1f5f9;
            color: var(--text-secondary);
        }}
        .badge-success {{ background: #d1fae5; color: #065f46; }}
        .badge-warning {{ background: #fef3c7; color: #92400e; }}
        .badge-danger {{ background: #fee2e2; color: #991b1b; }}

        /* ===== Thermo Section ===== */
        .thermo-grid {{
            display: grid;
            gap: 16px;
        }}
        .thermo-card {{
            background: #fafbfc;
            border: 1px solid var(--border);
            border-radius: var(--radius);
            overflow: hidden;
        }}
        .thermo-header {{
            background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
            padding: 16px 20px;
            display: flex;
            align-items: center;
            gap: 16px;
            flex-wrap: wrap;
            border-bottom: 1px solid var(--border);
        }}
        .thermo-rank {{
            background: var(--primary);
            color: white;
            padding: 4px 12px;
            border-radius: 6px;
            font-weight: 700;
            font-size: 12px;
        }}
        .thermo-seed {{
            font-family: monospace;
            font-size: 14px;
            color: var(--text);
        }}
        .thermo-seed code {{
            background: white;
            padding: 2px 8px;
            border-radius: 4px;
            border: 1px solid var(--border);
            font-weight: 600;
            color: var(--primary);
        }}
        .thermo-dg {{
            margin-left: auto;
            font-size: 14px;
            color: var(--text-secondary);
        }}
        .thermo-dg strong {{
            color: var(--text);
            font-weight: 700;
        }}
        .thermo-hits {{
            padding: 16px 20px;
        }}
        .thermo-subtitle {{
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--text-secondary);
            margin-bottom: 12px;
        }}
        .thermo-hits ul {{
            list-style: none;
        }}
        .thermo-hits li {{
            padding: 10px 14px;
            background: white;
            border-radius: 8px;
            margin-bottom: 8px;
            border: 1px solid var(--border);
            font-size: 13px;
        }}
        .hit-target {{
            font-weight: 600;
            color: var(--primary);
        }}
        .hit-meta {{
            color: var(--text-secondary);
            margin-left: 8px;
        }}
        .hit-struct {{
            display: block;
            margin-top: 6px;
            color: #94a3b8;
            font-family: monospace;
            background: #f8fafc;
            padding: 4px 8px;
            border-radius: 4px;
        }}

        .empty-state {{
            text-align: center;
            padding: 60px 20px;
            color: var(--text-secondary);
        }}
        .empty-icon {{
            font-size: 56px;
            margin-bottom: 16px;
        }}
        .empty-state h3 {{
            font-size: 20px;
            color: var(--text);
            margin-bottom: 8px;
        }}
        .empty-state p {{
            font-size: 15px;
            max-width: 500px;
            margin: 0 auto;
        }}

        /* ===== Parameters Grid ===== */
        .params-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 16px;
        }}
        .param-item {{
            display: flex;
            align-items: flex-start;
            gap: 12px;
            padding: 16px;
            background: #fafbfc;
            border-radius: 10px;
            border: 1px solid var(--border);
        }}
        .param-icon {{
            width: 36px;
            height: 36px;
            border-radius: 8px;
            background: #e0f2fe;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 16px;
            flex-shrink: 0;
        }}
        .param-content {{
            flex: 1;
        }}
        .param-label {{
            font-size: 12px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--text-secondary);
            margin-bottom: 4px;
        }}
        .param-value {{
            font-size: 15px;
            font-weight: 600;
            color: var(--text);
        }}
        .param-value code {{
            background: white;
            padding: 2px 6px;
            border-radius: 4px;
            border: 1px solid var(--border);
            font-family: monospace;
            font-size: 13px;
        }}

        /* ===== Footer ===== */
        .report-footer {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 32px 24px;
            border-top: 1px solid var(--border);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 12px;
            color: var(--text-secondary);
            font-size: 13px;
        }}
        .footer-brand {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-weight: 600;
            color: var(--text);
        }}
        .footer-links {{
            display: flex;
            gap: 20px;
        }}
        .footer-links a {{
            color: var(--text-secondary);
            text-decoration: none;
            transition: color 0.15s;
        }}
        .footer-links a:hover {{
            color: var(--primary);
        }}

        /* ===== Print ===== */
        @media print {{
            .report-header {{ padding: 24px; }}
            .section {{ break-inside: avoid; }}
            body {{ background: white; }}
            .stat-card, .section {{
                box-shadow: none;
                border: 1px solid #ddd;
            }}
        }}

        /* ===== Responsive ===== */
        @media (max-width: 768px) {{
            .report-header {{ padding: 32px 20px; }}
            .header-title {{ font-size: 26px; }}
            .stats-grid {{ grid-template-columns: 1fr 1fr; }}
            .section {{ padding: 20px; }}
            .thermo-header {{ flex-direction: column; align-items: flex-start; }}
            .thermo-dg {{ margin-left: 0; }}
        }}
    </style>
</head>
<body>
    <!-- Header -->
    <header class="report-header">
        <div class="header-content">
            <div class="header-badge">dsRNA-Forge Design Report</div>
            <h1 class="header-title">dsRNA-Forge 设计报告</h1>
            <p class="header-subtitle">昆虫/植物 dsRNA/siRNA/DsiRNA 与 SpCas9 sgRNA 设计工具 · 自动生成的实验验证报告</p>
            <div class="header-meta">
                <div class="header-meta-item">
                    <span>📅</span>
                    <span>生成时间: 2026-05-16</span>
                </div>
                <div class="header-meta-item">
                    <span>📦</span>
                    <span>版本: v0.1.1</span>
                </div>
                <div class="header-meta-item">
                    <span>⏱️</span>
                    <span>模式: {mode.upper()}</span>
                </div>
            </div>
        </div>
    </header>

    <!-- Stats Bar -->
    <div class="stats-bar">
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-icon mode">🧬</div>
                <div class="stat-value">{mode.upper()}</div>
                <div class="stat-label">设计模式</div>
            </div>
            <div class="stat-card">
                <div class="stat-icon total">📊</div>
                <div class="stat-value">{total}</div>
                <div class="stat-label">候选总数</div>
                <div class="stat-change success">Raw candidates generated</div>
            </div>
            <div class="stat-card">
                <div class="stat-icon passed">✅</div>
                <div class="stat-value">{passed}</div>
                <div class="stat-label">通过筛选</div>
                <div class="stat-change {'success' if passed > 0 else 'danger'}">{'All rules passed' if passed > 0 else 'No candidates passed'}</div>
            </div>
            <div class="stat-card">
                <div class="stat-icon vienna">🧪</div>
                <div class="stat-value">✓</div>
                <div class="stat-label">ViennaRNA</div>
                <div class="stat-change success">Thermodynamics enabled</div>
            </div>
        </div>
    </div>

    <!-- Main Content -->
    <main class="main-content">
        <!-- Results Table -->
        <section class="section">
            <div class="section-header">
                <div>
                    <h2 class="section-title">
                        <span class="section-title-icon">📋</span>
                        Top 候选结果
                    </h2>
                    <p class="section-subtitle">按 Consensus 评分排序的前 20 条候选序列，包含 GC 含量、脱靶风险与热力学评估</p>
                </div>
                <span class="section-badge">{min(total, 20)} / {total}</span>
            </div>
            <div class="table-wrapper">
                <table>
                    <thead>
                        <tr>
                            <th>Rank</th>
                            <th>Sequence</th>
                            <th>Consensus</th>
                            <th>GC%</th>
                            <th>Status</th>
                            <th>Risk</th>
                            <th>16bp</th>
                            <th>Seed7</th>
                            <th>ΔG (kcal/mol)</th>
                        </tr>
                    </thead>
                    <tbody>
                        {rows}
                    </tbody>
                </table>
            </div>
        </section>

        <!-- Thermodynamics -->
        <section class="section">
            <div class="section-header">
                <div>
                    <h2 class="section-title">
                        <span class="section-title-icon">🔥</span>
                        ViennaRNA 热力学脱靶分析
                    </h2>
                    <p class="section-subtitle">使用 RNAduplex 计算种子区 (nt 2-8) 与转录组的杂交自由能 · 阈值: -7.0 kcal/mol</p>
                </div>
            </div>
            <div class="thermo-grid">
                {thermo_details}
            </div>
        </section>

        <!-- Parameters -->
        <section class="section">
            <div class="section-header">
                <div>
                    <h2 class="section-title">
                        <span class="section-title-icon">⚙️</span>
                        设计参数
                    </h2>
                    <p class="section-subtitle">本次设计所使用的配置与规则</p>
                </div>
            </div>
            <div class="params-grid">
                <div class="param-item">
                    <div class="param-icon">🎯</div>
                    <div class="param-content">
                        <div class="param-label">目标基因</div>
                        <div class="param-value"><code>GeneA_CYP6G1_target</code></div>
                    </div>
                </div>
                <div class="param-item">
                    <div class="param-icon">📚</div>
                    <div class="param-content">
                        <div class="param-label">转录组</div>
                        <div class="param-value">3 个基因 (含 1 个潜在脱靶)</div>
                    </div>
                </div>
                <div class="param-item">
                    <div class="param-icon">📐</div>
                    <div class="param-content">
                        <div class="param-label">评分规则</div>
                        <div class="param-value">Consensus, Reynolds, Ui-Tei, Amarzguioui, Hsieh, Jagla</div>
                    </div>
                </div>
                <div class="param-item">
                    <div class="param-icon">🌡️</div>
                    <div class="param-content">
                        <div class="param-label">GC 范围</div>
                        <div class="param-value">30% — 55%</div>
                    </div>
                </div>
                <div class="param-item">
                    <div class="param-icon">🛡️</div>
                    <div class="param-content">
                        <div class="param-label">脱靶筛查</div>
                        <div class="param-value">16bp / 20bp / 7nt 种子区 + ViennaRNA 热力学</div>
                    </div>
                </div>
                <div class="param-item">
                    <div class="param-icon">🔧</div>
                    <div class="param-content">
                        <div class="param-label">工具链</div>
                        <div class="param-value">ViennaRNA RNAduplex / RNAcofold / RNAup</div>
                    </div>
                </div>
            </div>
        </section>
    </main>

    <!-- Footer -->
    <footer class="report-footer">
        <div class="footer-brand">
            <span>🧬</span>
            <span>dsRNA-Forge v0.1.1</span>
        </div>
        <div class="footer-links">
            <a href="https://github.com/dsRNA-Forge">GitHub</a>
            <a href="#">Documentation</a>
            <a href="#">Report Issue</a>
        </div>
    </footer>
</body>
</html>
"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    return output_path


if __name__ == "__main__":
    report_path = "/home/nee/dsRNA_Forge_MVP_review/dsRNA_Forge/report.html"
    generate_html(result, report_path)
    print(f"报告已生成: {report_path}")
