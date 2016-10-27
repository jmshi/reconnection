import IPython.core.display as di

def toggle_code():
  # This line will hide code by default when the notebook is exported as HTML
  di.display_html('<script>jQuery(function() {if (jQuery("body.notebook_app").length == 0) { jQuery(".input_area").toggle(); jQuery(".prompt").toggle();}});</script>', raw=True)
  
  # This line will add a button to toggle visibility of code blocks, for use with the HTML export version
  di.display_html('''<button onclick="jQuery('.input_area').toggle(); jQuery('.prompt').toggle();">Toggle code</button>''', raw=True)

  
if __name__=='__main__':
    if len( sys.argv ) > 0:
        print "Please do not specify input"
        exit(  );
    #

    toggle_code();
    #shwave.vx_amp('../develop/bin/debug',0,400)
#

