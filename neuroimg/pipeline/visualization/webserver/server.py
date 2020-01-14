from flask import Flask, render_template, Response, send_file
import os


app = Flask(__name__)

app._static_folder = os.path.abspath("templates/static")


@app.route("/")
def home():
    return render_template("threeD.html")

  
@app.route('/brain')
def api_articles():
	  return send_file('./templates/static/reconstruction.glb')

if __name__ == "__main__":
    app.run(debug=True)

