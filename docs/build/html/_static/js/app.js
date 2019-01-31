jQuery(function($){
	var $docsmenu = $('aside > ul');
	
	// hide all subs first
	showcurrent = function (){
	    
	    $docsmenu.find('li.current > ul').show();
	    
	    var anchor = encodeURI(window.location.hash);
	    if (!anchor) return;
	    try {
		var link = $('aside').find('[href="' + anchor + '"]');

		link.addClass('current');
		
		link.closest('li.toctree-l3').addClass('current');
		link.closest('li.toctree-l3 ul').addClass('current');
		link.closest('li.toctree-l4').addClass('current');

	    }
	    catch (err) {
		console.log('Error expanding nav for anchor', err);
	    }

	}
	showcurrent();
	
	$docsmenu.find('.toctree-l3:not(.current) ul,.toctree-l2:not(.current) ul,.toctree-l1:not(.current) ul').hide();
		
	$docsmenu.find('a').click(function() {
	    if (!/^#/.test($(this).attr('href'))) return;
	    
	    $(this).closest('ul').find('li').removeClass('current');
	    $(this).closest('ul').find('a').removeClass('current');
	    $(this).addClass('current');
	    $(this).parent().addClass('current');
	});
	
	// loop through clickable items
	$docsmenu.find('.toctree-l3 a,.current .toctree-l2 a').each(function() { 
		$(this).click(function() {
		
			// slide subs down and up
			$(this).next('ul').slideToggle(200);
			$(this).next('ul').siblings('a').toggleClass('expanded');
		
			// hide other subs when new one is clicked
			$(this).parent().siblings().children().next().slideUp();
			$(this).parent().siblings().children().removeClass('expanded');
		});


		// if item has subs, give class
		if ($(this).siblings().size() > 0) { 
			$(this).addClass('has-subs');
			
			// if item has subs, add toggle button
			$(this).prepend('<span class="toggle-btn"></span>');
			if( $(this).next('ul').hasClass('current') ){
			    $(this).next('ul').siblings('a').toggleClass('expanded');
			}
		} 

		// if item has no-subs, give different class
		else { 
			$(this).addClass('has-no-subs');
		}


		// prevent click from jumping to link
		$('.toggle-btn').click(function(event) {
			event.preventDefault();
		});
	}
    );
});


jQuery(function($){
	// get hash from url
	var hash = window.location.hash;

	if(hash.length != 0) {
		// expands accordeon panel when page is entered through anchorlink
		var idToToggle = hash.replace('#', '');
		$('#'+idToToggle).collapse('show');
	}
});


jQuery(function($){
	$('input#search-faq').quicksearch('.faq-item', {
		delay: 180,
		onBefore: function(){
			$('.panel-collapse.in').removeClass('in');
			$('.hilite1').contents().unwrap();
		},
		loaderText: 'Searching...',
		onAfter: function(){
			if($(this).val() != '') {
				$('.faq-item:visible .panel-collapse').addClass('in');
			}

			$('.faq-item:visible .panel-collapse .answer').SearchHighlight({ exact: 'exact', highlight: $('.faq-item:visible .panel-collapse .answer'), keys: $('input#search-faq').val() });
		}
	});
});

var scm_hideFooter = false;

function showHideFooter(){

	if (jQuery(window).width() < 430) {
		// remove class expanded on all .widget_nav_menu h3
		jQuery('.widget_nav_menu h3').removeClass('expanded');

		// add to all class hidden-xs
		jQuery('.widget_nav_menu div').addClass('hidden-xs');

		jQuery('.widget_nav_menu h3').click(function() {

				// this
				if( jQuery(this).hasClass('expanded') ) {
					var isExpanded = true;
				} else {
					var isExpanded = false;
				}

				// remove class expanded on all .widget_nav_menu h3
				jQuery('.widget_nav_menu h3').removeClass('expanded');

				// add to all class hidden-xs
				jQuery('.widget_nav_menu div').addClass('hidden-xs');

				if( !isExpanded ) {

					// add expanded .widget_nav_menu h3
					jQuery(this).addClass('expanded');
					// remove class hidden-xs
					jQuery(this).next('div').removeClass('hidden-xs');

				}
		} );
	} else {

		// remove class expanded on all .widget_nav_menu h3
		jQuery('.widget_nav_menu h3').removeClass('expanded');

		// add to all class hidden-xs
		jQuery('.widget_nav_menu div').removeClass('hidden-xs');
	}
}

// load on resize and on load
jQuery(window).resize(function() {
	showHideFooter();
});

jQuery(function(){
	showHideFooter();
});


