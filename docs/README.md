# GELATOドキュメント

このディレクトリには、GELATOプロジェクトのSphinxドキュメントが含まれています。

## ドキュメントのビルド

### 必要なパッケージのインストール

```bash
pip install -r requirements.txt
```

### HTMLドキュメントのビルド

```bash
make html
```

生成されたHTMLドキュメントは `_build/html/` ディレクトリに出力されます。

### ブラウザで表示

```bash
# Linux/macOS
open _build/html/index.html

# Windows
start _build/html/index.html
```

## その他のフォーマット

Sphinxは様々な出力フォーマットをサポートしています：

- `make html` - HTMLドキュメント
- `make latexpdf` - PDFドキュメント（LaTeXが必要）
- `make epub` - ePub形式
- `make man` - manページ
- `make text` - プレーンテキスト

全ての利用可能なターゲットを確認するには：

```bash
make help
```

## ドキュメント構成

- `conf.py` - Sphinx設定ファイル
- `index.rst` - トップページ
- `installation.rst` - インストール手順
- `usage.rst` - 使用方法
- `modules.rst` - モジュール一覧
- `api.rst` - API完全リファレンス
- `api/` - 個別モジュールのAPIドキュメント

## 自動ドキュメント生成

このドキュメントは `sphinx.ext.autodoc` を使用してPythonのdocstringから自動生成されます。
コード内にdocstringを追加・編集した場合は、再ビルドすることでドキュメントに反映されます。
