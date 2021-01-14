library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity full_regionizer_mux_stream is
    generic(
            MU_ETA_CENTER : integer 
    );
    port(
            ap_clk : IN STD_LOGIC;
            ap_rst : IN STD_LOGIC;
            ap_start : IN STD_LOGIC;
            ap_done : OUT STD_LOGIC;
            ap_idle : OUT STD_LOGIC;
            ap_ready : OUT STD_LOGIC;
            tracks_start    : IN STD_LOGIC;
            tracks_newevent : IN STD_LOGIC;
            tracks_in_0_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_0_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_1_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_1_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_2_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_2_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_3_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_3_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_4_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_4_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_5_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_5_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_6_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_6_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_7_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_7_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_8_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_8_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_start    : IN STD_LOGIC;
            calo_newevent : IN STD_LOGIC;
            calo_in_0_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            mu_start    : IN STD_LOGIC;
            mu_newevent : IN STD_LOGIC;
            mu_in_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            mu_in_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            --
            tracks_out : OUT w64s(NTKSTREAM-1 downto 0);
            calo_out   : OUT w64s(NCALOSTREAM-1 downto 0);
            mu_out     : OUT w64s(NMUSTREAM-1   downto 0);
            newevent_out : OUT STD_LOGIC
    );
end full_regionizer_mux_stream;

architecture Behavioral of full_regionizer_mux_stream is

    signal tracks_regionized:        w64s(NPFREGIONS-1 downto 0);
    signal tracks_regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal tracks_regionized_roll:   std_logic := '0';

    signal calo_regionized:        w64s(NPFREGIONS-1 downto 0);
    signal calo_regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal calo_regionized_roll:   std_logic := '0';

    signal mu_regionized:        w64s(NPFREGIONS-1 downto 0);
    signal mu_regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal mu_regionized_roll:   std_logic := '0';

begin

    tk_regionizer : entity work.tk_regionizer 
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => tracks_start,
                             newevent => tracks_newevent,
                             tracks_in_0_0_V => tracks_in_0_0_V,
                             tracks_in_0_1_V => tracks_in_0_1_V,
                             tracks_in_1_0_V => tracks_in_1_0_V,
                             tracks_in_1_1_V => tracks_in_1_1_V,
                             tracks_in_2_0_V => tracks_in_2_0_V,
                             tracks_in_2_1_V => tracks_in_2_1_V,
                             tracks_in_3_0_V => tracks_in_3_0_V,
                             tracks_in_3_1_V => tracks_in_3_1_V,
                             tracks_in_4_0_V => tracks_in_4_0_V,
                             tracks_in_4_1_V => tracks_in_4_1_V,
                             tracks_in_5_0_V => tracks_in_5_0_V,
                             tracks_in_5_1_V => tracks_in_5_1_V,
                             tracks_in_6_0_V => tracks_in_6_0_V,
                             tracks_in_6_1_V => tracks_in_6_1_V,
                             tracks_in_7_0_V => tracks_in_7_0_V,
                             tracks_in_7_1_V => tracks_in_7_1_V,
                             tracks_in_8_0_V => tracks_in_8_0_V,
                             tracks_in_8_1_V => tracks_in_8_1_V,
                             tracks_out_0_V => tracks_regionized(0),
                             tracks_out_1_V => tracks_regionized(1),
                             tracks_out_2_V => tracks_regionized(2),
                             tracks_out_3_V => tracks_regionized(3),
                             tracks_out_4_V => tracks_regionized(4),
                             tracks_out_5_V => tracks_regionized(5),
                             tracks_out_6_V => tracks_regionized(6),
                             tracks_out_7_V => tracks_regionized(7),
                             tracks_out_8_V => tracks_regionized(8),
                             tracks_out_valid_0 => tracks_regionized_valid(0),
                             tracks_out_valid_1 => tracks_regionized_valid(1),
                             tracks_out_valid_2 => tracks_regionized_valid(2),
                             tracks_out_valid_3 => tracks_regionized_valid(3),
                             tracks_out_valid_4 => tracks_regionized_valid(4),
                             tracks_out_valid_5 => tracks_regionized_valid(5),
                             tracks_out_valid_6 => tracks_regionized_valid(6),
                             tracks_out_valid_7 => tracks_regionized_valid(7),
                             tracks_out_valid_8 => tracks_regionized_valid(8),
                             newevent_out => tracks_regionized_roll);


    calo_regionizer : entity work.calo_regionizer 
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => calo_start,
                             newevent => calo_newevent,
                             calo_in_0_0_V => calo_in_0_0_V,
                             calo_in_0_1_V => calo_in_0_1_V,
                             calo_in_0_2_V => calo_in_0_2_V,
                             calo_in_0_3_V => calo_in_0_3_V,
                             calo_in_1_0_V => calo_in_1_0_V,
                             calo_in_1_1_V => calo_in_1_1_V,
                             calo_in_1_2_V => calo_in_1_2_V,
                             calo_in_1_3_V => calo_in_1_3_V,
                             calo_in_2_0_V => calo_in_2_0_V,
                             calo_in_2_1_V => calo_in_2_1_V,
                             calo_in_2_2_V => calo_in_2_2_V,
                             calo_in_2_3_V => calo_in_2_3_V,
                             calo_out_0_V => calo_regionized(0),
                             calo_out_1_V => calo_regionized(1),
                             calo_out_2_V => calo_regionized(2),
                             calo_out_3_V => calo_regionized(3),
                             calo_out_4_V => calo_regionized(4),
                             calo_out_5_V => calo_regionized(5),
                             calo_out_6_V => calo_regionized(6),
                             calo_out_7_V => calo_regionized(7),
                             calo_out_8_V => calo_regionized(8),
                             calo_out_valid_0 => calo_regionized_valid(0),
                             calo_out_valid_1 => calo_regionized_valid(1),
                             calo_out_valid_2 => calo_regionized_valid(2),
                             calo_out_valid_3 => calo_regionized_valid(3),
                             calo_out_valid_4 => calo_regionized_valid(4),
                             calo_out_valid_5 => calo_regionized_valid(5),
                             calo_out_valid_6 => calo_regionized_valid(6),
                             calo_out_valid_7 => calo_regionized_valid(7),
                             calo_out_valid_8 => calo_regionized_valid(8),
                             newevent_out => calo_regionized_roll);

    mu_regionizer : entity work.mu_regionizer 
                generic map(ETA_CENTER => MU_ETA_CENTER)
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => mu_start,
                             newevent => mu_newevent,
                             mu_in_0_V => mu_in_0_V,
                             mu_in_1_V => mu_in_1_V,
                             mu_out_0_V => mu_regionized(0),
                             mu_out_1_V => mu_regionized(1),
                             mu_out_2_V => mu_regionized(2),
                             mu_out_3_V => mu_regionized(3),
                             mu_out_4_V => mu_regionized(4),
                             mu_out_5_V => mu_regionized(5),
                             mu_out_6_V => mu_regionized(6),
                             mu_out_7_V => mu_regionized(7),
                             mu_out_8_V => mu_regionized(8),
                             mu_out_valid_0 => mu_regionized_valid(0),
                             mu_out_valid_1 => mu_regionized_valid(1),
                             mu_out_valid_2 => mu_regionized_valid(2),
                             mu_out_valid_3 => mu_regionized_valid(3),
                             mu_out_valid_4 => mu_regionized_valid(4),
                             mu_out_valid_5 => mu_regionized_valid(5),
                             mu_out_valid_6 => mu_regionized_valid(6),
                             mu_out_valid_7 => mu_regionized_valid(7),
                             mu_out_valid_8 => mu_regionized_valid(8),
                             newevent_out => mu_regionized_roll);

    tk_delay_sort_mux_stream : entity work.delay_sort_mux_stream
                generic map(NREGIONS => NPFREGIONS, 
                            NSORTED  => NTKSORTED,
                            NSTREAM  => NTKSTREAM,
                            OUTII    => PFII240,
                            DELAY    => TKDELAY)
                port map(ap_clk => ap_clk,
                         d_in => tracks_regionized,
                         valid_in => tracks_regionized_valid,
                         roll => tracks_regionized_roll,
                         d_out => tracks_out,
                         roll_out => newevent_out);

    calo_delay_sort_mux_stream : entity work.delay_sort_mux_stream
                generic map(NREGIONS => NPFREGIONS, 
                            NSORTED  => NCALOSORTED,
                            NSTREAM  => NCALOSTREAM,
                            OUTII    => PFII240,
                            DELAY    => CALODELAY)
                port map(ap_clk => ap_clk,
                         d_in => calo_regionized,
                         valid_in => calo_regionized_valid,
                         roll => calo_regionized_roll,
                         d_out => calo_out,
                         roll_out => open);

    mu_delay_sort_mux_stream : entity work.delay_sort_mux_stream
                generic map(NREGIONS => NPFREGIONS, 
                            NSORTED  => NMUSORTED,
                            NSTREAM  => NMUSTREAM,
                            OUTII    => PFII240,
                            DELAY    => MUDELAY)
                port map(ap_clk => ap_clk,
                         d_in => mu_regionized,
                         valid_in => mu_regionized_valid,
                         roll => mu_regionized_roll,
                         d_out => mu_out,
                         roll_out => open);
end Behavioral;
